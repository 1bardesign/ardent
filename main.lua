local ffi = require("ffi")
require("lib.core")
require("matrix")
lg = love.graphics
love.window.setFullscreen(true)

local vert_res = 64
local vert_count = vert_res * vert_res

local verts = {}
for y = 1, vert_res do
	for x = 1, vert_res do
		table.insert(verts, {
			(x - 0.5) / vert_res,
			(y - 0.5) / vert_res,
		})
	end
end
local vert_mesh = lg.newMesh(
	{{"VertexUV", "float", 2}},
	verts,
	"points",
	"static"
)
local vert_mesh_shader = lg.newShader([[
uniform Image position_buffer;
uniform Image colour_map;
uniform bool draw_colour;

uniform vec2 res;

uniform mat4 view;
uniform mat4 proj;

varying float v_col;
varying float v_depth;

#ifdef VERTEX
attribute vec2 VertexUV;
vec4 position(mat4 love_proj, vec4 _pos) {
	vec4 t = Texel(position_buffer, VertexUV);
	vec3 pos = t.rgb;
	v_col = t.a;

	vec4 r = proj * view * vec4(pos, 1.0);

	r.x /= res.x / res.y;

	v_depth = 1.0;// - r.z;

	return r;
}
#endif
#ifdef PIXEL
vec4 effect(vec4 col, Image tex, vec2 uv, vec2 screenpos) {
	if (!draw_colour) {
		return vec4(1.0);
	}
	return Texel(colour_map, vec2(v_col, 0.0)) * v_depth;
}
#endif
]])

function init_view()
	view = {
		dirty = 	true,
		from = 		vec3:xyz(0, 0, -1):smuli(50),
		to = 		vec3:zero(),
		right =     vec3:xyz(1,0,0),
		up =        vec3:xyz(0,-1,0),
		forward =   vec3:xyz(0,0,1),
		zoom = 		vec3:xyz(1,1,1):smuli(0.01),
	}
end

function update_cam(dt)

	if not _g_old_view then
		_g_old_view = {}
		--copy the vectors
		for k,v in pairs(view) do
			if type(v) == "table" and v.copy then
				_g_old_view[k] = v:copy()
			end
		end
	end

	--zoom
	local zoom_factor = 1 + dt
	if love.keyboard.isDown("z") then view.zoom:smuli(zoom_factor) end
	if love.keyboard.isDown("x") then view.zoom:sdivi(zoom_factor) end

	--rotate
	local turn_speed = math.pi * 0.5 * dt
	local function rotate_axes(on_axis, angle)
		local m = rotation_axis_angle(on_axis, angle)
		for i,v in ipairs {
			"right", "up", "forward", "from",
		} do
			view[v] = transform_direction(view[v], m)
		end
	end
	if love.keyboard.isDown("w") then rotate_axes(view.right, turn_speed) end
	if love.keyboard.isDown("s") then rotate_axes(view.right, -turn_speed) end
	if love.keyboard.isDown("a") then rotate_axes(view.up, turn_speed) end
	if love.keyboard.isDown("d") then rotate_axes(view.up, -turn_speed) end
	if love.keyboard.isDown("q") then rotate_axes(view.forward, turn_speed) end
	if love.keyboard.isDown("e") then rotate_axes(view.forward, -turn_speed) end

	--todo: re-orthoganalise

	--check for changes automagically
	local changed = false
	for k,v in pairs(_g_old_view) do
		if v:nequals(view[k]) then
			changed = true
			break
		end
	end
	if changed then
		view.dirty = true
		for k,v in pairs(_g_old_view) do
			v:vset(view[k])
		end
	end
	if view.dirty then
		
		--translation
		local sx, sy, sz = view.zoom:unpack()
		local scale_mat = {
			sx,0,0,0,
			0,sy,0,0,
			0,0,sz,0,
			0,0,0,1,
		}
		local cam_mat = calculate_camera(view.from, view.to, view.up)
		local view_mat = matmul4x4(cam_mat, scale_mat)

		vert_mesh_shader:send("res", {lg:getDimensions()})
		vert_mesh_shader:send("view", view_mat)
		-- vert_mesh_shader:send("proj", calculate_perspective({lg:getDimensions()}, 1.0, 0.01, 200))
		vert_mesh_shader:send("proj", {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1,})
	end
end

local transform_functions = {
	{"linear", ([[
		t = i;
	]]),},
	--multi component systems
	{"sin", ([[
		t = sin(i);
	]]),},
	--single component waves
	{"sinx", ([[
		t = vec3(i.x, sin(i.x), i.z);
	]]),},
	--sqrt all components
	{"sqrt", ([[
		t = sign(i) * sqrt(abs(i));
	]]),},
	--classical flame transforms extended to 3d where possible
	{"spherical", ([[
		t = i / r2;
	]]),},
	{"swirl", ([[
		float c = cos(r2);
		float s = sin(r2);
		t = vec3(
			c * i.x - s * i.y,
			s * i.x + c * i.y,
			i.z
		);
	]]),},
	{"horseshoe", ([[
		t = vec3(
			(i.x - i.y) * (i.x + i.y),
			2.0 * i.x * i.y,
			i.z
		) / r;
	]]),},
	{"polar", ([[
		t = vec3(
			theta / pi,
			r - 1.0,
			i.z
		);
	]])},
}

--autogen cdef
ffi.cdef([[
struct flame_transform{
	//transforms
	float pre[12];
	float post[12];

	//function lerps
	float function_weights[]]..#transform_functions..[[];
};
]])

local transform_functions_ordered = {}
for i,v in ipairs(transform_functions) do
	table.insert(transform_functions_ordered, v[1])
end

local functions_snippet = {}
for i,v in ipairs(transform_functions) do
	local var_name = "f_"..v[1]
	table.insert(functions_snippet, table.concat{
		"float ", var_name, " = transforms[idx++];\n",
		"if (abs(", var_name, ") != 0.0) {\n",
		"\tvec3 t = vec3(0.0);",
		v[2], "\n",
		"\to = ", var_name, " * t + o;\n",
		"}",
	})
end
functions_snippet = table.concat(functions_snippet, "\n")

local position_buffer = lg.newCanvas(vert_res, vert_res, {format="rgba32f"})
local position_buffer_previous = lg.newCanvas(vert_res, vert_res, {format="rgba32f"})
local position_shader = lg.newShader([[
const int MAX_TRANSFORMS = 8;
const int TRANSFORM_SIZE = 12 + 12 + 6;
uniform int transform_count;
uniform int transform_offset;
uniform int transform_stride;
uniform float transforms[MAX_TRANSFORMS * TRANSFORM_SIZE];

mat4 read_mat4(int idx) {
	return mat4(
		transforms[idx+0], transforms[idx+4], transforms[idx+8],  0.0,
		transforms[idx+1], transforms[idx+5], transforms[idx+9],  0.0,
		transforms[idx+2], transforms[idx+6], transforms[idx+10], 0.0,
		transforms[idx+3], transforms[idx+7], transforms[idx+11], 1.0
	);
}

vec3 tr(vec3 i, int tr) {
	int idx = tr * TRANSFORM_SIZE;

	///////////////////////////////////
	//collect transforms
	mat4 a = read_mat4(idx); idx += 12;
	mat4 b = read_mat4(idx); idx += 12;

	///////////////////////////////////
	//pre-transform
	i = (a * vec4(i, 1.0)).xyz;

	///////////////////////////////////
	//apply activation

	//common parameters
	//distance from origin
	float r2 = i.x * i.x + i.y * i.y + i.z * i.z;
	float r = sqrt(r2);
	
	//angles; todo 3d?
	float theta = atan(i.y / i.x);
	float phi = atan(i.x / i.y);
	const float pi = 3.14159265359;

	//accumulate output in here
	vec3 o = vec3(0.0);
	
	]]..functions_snippet..[[

	///////////////////////////////////
	//post-transform
	o = (b * vec4(o, 1.0)).xyz;

	return o;
}

#ifdef PIXEL
vec4 effect(vec4 color, Image tex, vec2 uv, vec2 screenpos) {
	//pick the transform
	float transform = floor(mod(
		floor(
			float(int(screenpos.y) / transform_stride)
			+ float(transform_offset)
		),
		float(transform_count)
	));
	//sample the point
	vec4 t = Texel(tex, uv);
	//apply
	return vec4(
		tr(t.xyz, int(transform)),
		//equidistant colours
		mix(t.w, float(transform) / (float(transform_count) - 1.0), 0.5)
	);
}
#endif
]])

local transforms = {}

function to_tr(t)
	--todo: parse nice transforms again
	if true then
		return nil
	end

	local c = math.cos(t.angle)
	local s = math.sin(t.angle)

	local ox, oy = unpack(t.offset)
	local ax, ay = unpack(t.around)
	local sx, sy = unpack(t.scale)

	local rx = ax - c * ax + s * ay
	local ry = ay - s * ax - c * ay

	return {
		c * sx, -s,     ox + rx,
		s,      c * sy, oy + ry,
	}
end

function upload_transforms()
	local count = #transforms
	local tr_d = love.data.newByteData(ffi.sizeof("struct flame_transform") * count)
	local tr_ser = ffi.cast("struct flame_transform*", tr_d:getFFIPointer())
	for i,v in ipairs(transforms) do
		local t = tr_ser[i - 1]

		--copy rule amounts
		local rule_amounts = {}
		local func = v.func
		for i, v in ipairs(transform_functions_ordered) do
			local f = 0
			if type(func) == "table" then
				if func[v] then
					--key value table
					f = func[v]
				else
					--ordered table of either
					--	{"f", v}
					--	or "f"
					--	or in-order numeric rule weights
					for ci, cf in ipairs(func) do
						if type(cf) == "string" and (v == cf) then
							f = 1
						elseif type(cf) == "table" and v == cf[1] then
							f = v[2]
						elseif type(cf) == "number" and ci == i then
							f = cf
						end
					end
				end
			elseif type(func) == "string" then
				f = (func == v) and 1 or 0
			elseif type(func) == "number" then
				f = (func == i) and 1 or 0
			end
			table.insert(rule_amounts, f)
		end

		--copy transforms
		for i = 1, 12 do
			t.pre[i  - 1] = v.pre[i]
			t.post[i - 1] = v.post[i]
		end

		--todo: this could probably be an array :)
		for i,v in ipairs(rule_amounts) do
			t.function_weights[i - 1] = v
		end
	end
	position_shader:send("transforms", tr_d)
	position_shader:send("transform_count", count)
end

function random_transforms()
	--
	transforms = {}

	local function _random_rule()
		return love.math.random(1, #transform_functions)
	end

	--primary transform for the set
	local primary_transform = _random_rule()

	--generate a certain number of rules
	for i = 1, love.math.random(3, 6) do
		--randomly select a few rules with random weighting
		local rules_max = love.math.random(1, 3)
		local rules_count  = love.math.random(1, rules_max)
		local function_weights = {}
		for i,v in ipairs(transform_functions_ordered) do function_weights[i] = 0 end
		for i = 0, rules_count do
			local f = _random_rule()
			if i == 0 then
				f = primary_transform
			end
			function_weights[f] = love.math.randomNormal()
			print(transform_functions[f][1], function_weights[f])
		end
		--normalise
		local t = 0
		for i,v in ipairs(function_weights) do
			t = t + math.abs(function_weights[i])
		end
		for i,v in ipairs(function_weights) do
			function_weights[i] = v / t
		end

		local function compute_transform(tr)
			--start with identity
			local t = {
				1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1,
			}
			--scale
			local sx, sy, sz = unpack(tr.scale)
			local scale_mat = {
				sx, 0, 0, 0,
				0, sy, 0, 0,
				0, 0, sz, 0,
				0, 0, 0, 1,
			};

			--rotate
			local rotation_mat = {
				1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1,
			}
			if t.euler then
				--euler angles
				do
					local a = t.euler[1]
					local s = math.sin(a)
					local c = math.cos(a)
					rotation_mat = matmul4x4(rotation_mat, {
						c, -s, 0, 0,
						s, c, 0, 0,
						0, 0, 1, 0,
						0, 0, 0, 1,
					})
				end
				do
					local a = t.euler[1]
					local s = math.sin(a)
					local c = math.cos(a)
					rotation_mat = matmul4x4(rotation_mat, {
						c, 0, s, 0,
						0, 1, 0, 0,
						-s, 0, c, 0,
						0, 0, 0, 1,
					})
				end
				do
					local a = t.euler[3]
					local s = math.sin(a)
					local c = math.cos(a)
					rotation_mat = matmul4x4(rotation_mat, {
						1, 0, 0, 0,
						0, c, -s, 0,
						0, s, c, 0,
						0, 0, 0, 1,
					})
				end
			elseif t.quat then
				--quaternion
				local w, x, y, z = unpack(t.quat)

				--cayley
				if w ~= 0 then
					x = x / w
					y = y / w
					z = z / w
					w = w / w
				end

				local n = w * w + x * x + y * y + z * z

				--normalize
				if false and n ~= 0 then
					local len = math.sqrt(n)
					w = w / len
					x = x / len
					y = y / len
					z = z / len

					n = 1
				end


				local s = 0
				if n ~= 0 then
					--
					s = 2 / n

				end
				local wx = s * w * x
				local wy = s * w * y
				local wz = s * w * z

				local xx = s * x * x
				local xy = s * x * y
				local xz = s * x * z

				local yy = s * y * y
				local yz = s * y * z

				local zz = s * z * z

				rotation_mat = {
					1 - (yy + zz), 	xy - wz,		xz + wy, 		0,
					xy + wz, 		1 - (xx + zz),	yz - wx,		0,
					xz - wy,		yz + wx,		1 - (xx + yy),	0,
					0,				0,				0,				1,
				}

			end

			--apply offset
			local ox, oy, oz = unpack(tr.offset)
			local offset_mat = {
				1, 0, 0, ox,
				0, 1, 0, oy,
				0, 0, 1, oz,
				0, 0, 0, 1,
			}

			t = matmul4x4(scale_mat, t)
			t = matmul4x4(rotation_mat, t)
			t = matmul4x4(offset_mat, t)


			--remove extra line which we don't upload
			for i = 1, 4 do
				table.remove(t)
			end

			return t
		end

		local function rt()
			local function r(o, s)
				return love.math.randomNormal(o, s)
			end
			return {
				offset = {
					r(0, 1), r(0, 1), r(0, 1),
				},
				scale = {
					r(1, 0.25), r(1, 0.25), r(1, 0.25),
				},
				-- euler = {
				-- 	r(0, math.pi), r(0, math.pi), r(0, math.pi),
				-- }
				quat = {
					r(0, 1), r(0, 1), r(0, 1), r(0, 1),
				},
			}
		end

		local t = {
			pre = compute_transform(rt()),
			post = compute_transform(rt()),
			func = function_weights,
		}
		table.insert(transforms, t)
	end
	upload_transforms()
end

random_transforms()

function random_positions()
	position_buffer:setFilter("nearest", "nearest")
	position_buffer_previous:setFilter("nearest", "nearest")

	local id = love.image.newImageData(vert_res, vert_res, "rgba32f")
	id:mapPixel(function(x, y, r, g, b, a)
		return
			--pos
			love.math.random() * 2 - 1,
			love.math.random() * 2 - 1,
			love.math.random() * 2 - 1,
			--col
			love.math.random()
	end)
	local im = love.graphics.newImage(id)
	lg.setCanvas(position_buffer)
	lg.setBlendMode("replace", "premultiplied")
	lg.draw(im)
	lg.setBlendMode("alpha")
	lg.setCanvas()

	upload_transforms()
end

function update_positions()
	--swap buffers
	position_buffer, position_buffer_previous = position_buffer_previous, position_buffer
	lg.setBlendMode("replace", "premultiplied")
	lg.setCanvas(position_buffer)
	lg.setShader(position_shader)
	--update random transform offset per-row
	local offset = love.math.random(0, #transforms)
	position_shader:send("transform_offset", offset)
	position_shader:send("transform_stride", love.math.random(1, 5))

	--rotated to get a bit more scrambling across elements
	lg.draw(
		position_buffer_previous,
		vert_res, 0,
		math.pi * 0.5
	)
	lg.setCanvas()
	lg.setShader()
	lg.setBlendMode("alpha")
end

local max_downres_factor = 4
local max_shader = lg.newShader([[
const int sample_size = ]]..tostring(max_downres_factor)..[[;
uniform vec2 res;
#ifdef PIXEL
vec4 effect(vec4 color, Image tex, vec2 uv, vec2 screenpos) {
	//maximum
	vec4 maximum = vec4(0.0);
	for (int y = 0; y < sample_size; y++)
	{
		for (int x = 0; x < sample_size; x++)
		{
			vec2 _uv = 
				//base uv
				floor(uv * res)
				//offset to edge
				- vec2(sample_size * 0.5)
				//offset per-pixel
				+ vec2(x, y)
				//center of texel
				+ vec2(0.5);
			maximum = max(maximum, Texel(tex, _uv / res));
		}
	}
	return maximum;
}
#endif
]])
function calc_max(t)
	local w, h = t:getDimensions()
	local fmt_t = {format=t:getFormat()}
	local e_fac = 1 / max_downres_factor
	local junk_cv = {}
	lg.setBlendMode("none")
	lg.setShader(max_shader)
	while w > 1 and h > 1 do
		local nw = math.ceil(w * e_fac)
		local nh = math.ceil(h * e_fac)
		local nt = love.graphics.newCanvas(nw, nh, fmt_t)
		table.insert(junk_cv, nt)
		lg.setCanvas(nt)
		max_shader:send("res", {t:getDimensions()})
		t:setFilter("nearest", "nearest")
		lg.draw(t, 0, 0, 0, e_fac, e_fac)
		lg.setCanvas()
		t = nt
		w, h = t:getDimensions()
	end
	lg.setBlendMode("alpha")
	lg.setShader()
	lg.setCanvas()
	--clean up sooner than later
	local final = table.remove(junk_cv)
	for i,v in ipairs(junk_cv) do
		v:release()
	end
	return final
end


local sw, sh = lg.getDimensions()
--todo: supersample
local screen_buffer = lg.newCanvas(
	sw, sh,
	{format="rgba32f"}
)
--uuuurgh we cant add on alpha channel currently
local count_buffer = lg.newCanvas(
	sw, sh,
	{format="r32f"}
)

local screen_shader = lg.newShader([[
uniform Image maximum;
uniform Image count;
#ifdef PIXEL
vec4 effect(vec4 color, Image tex, vec2 uv, vec2 screenpos) {
	float highest = Texel(maximum, vec2(0.5)).a;

	vec4 t = Texel(tex, uv);
	//replace alpha with count buffer
	t.a = Texel(count, uv).r;

	//mean
	vec3 rgb = t.rgb / t.a;
	//linear
	//float alpha = t.a / highest;
	//loglog
	float alpha = log(1.0 + t.a) / log(1.0 + highest);

	//power / gamma
	alpha = pow(alpha, 1.0 / 2.2);

	//scale
	alpha *= 0.5;

	return vec4(
		rgb * alpha,
		1.0
	);
}
#endif
]])

function hsl_to_rgb(h, s, l)
	local c = (1 - math.abs(2 * l - 1)) * s
	local hp = (h % 1) * 6
	local x = (1 - math.abs(hp % 2 - 1)) * c
	local r, g, b = 0, 0, 0
	if hp < 1 then
		r, g, b = c, x, 0
	elseif hp < 2 then
		r, g, b = x, c, 0
	elseif hp < 3 then
		r, g, b = 0, c, x
	elseif hp < 4 then
		r, g, b = 0, x, c
	elseif hp < 5 then
		r, g, b = x, 0, c
	elseif hp < 6 then
		r, g, b = c, 0, x
	end
	local m = l - c / 2
	r = r + m
	g = g + m
	b = b + m
	return r, g, b
end

local colour_map
function new_colours()
	local res = 256
	local id = love.image.newImageData(res, 1)
	
	--hue params
	local h_o = love.math.random()
	local h_g = love.math.randomNormal(0, 0.1)

	--saturation oscilator
	local s_o = love.math.random() * math.pi * 2
	local s_f = love.math.random() * math.pi * 2
	local s_b = love.math.random() * 0.5 + 0.25
	local s_a = love.math.random()

	--lightness fbm
	local l = love.math.random()
	local l_jump_chance = 10 / res
	local l_wander = 20 / res

	--noise
	local n_mono = love.math.random() * 0.1
	local n_chro = love.math.random() * 0.1

	local min_l = 0.1
	local max_l = 0.5

	id:mapPixel(function(x, y, r, g, b, a)
		local f = (x / res) * 2 * math.pi

		local h = h_o + f * h_g
		local s = s_b + s_a * math.sin(s_o + s_f * f) * 0.5
		--wander lightness
		if love.math.random() < l_jump_chance then
			l = (l + love.math.random() * 2 - 1)
		end
		l = l + (love.math.random() * 2 - 1) * l_wander
		l = math.max(min_l, math.min(l, max_l))

		local r, g, b = hsl_to_rgb(h, s, l)

		--add noise
		local mono = (love.math.random() - 0.5) * n_mono
		r = r + ((love.math.random() - 0.5) * n_chro + mono)
		g = g + ((love.math.random() - 0.5) * n_chro + mono)
		b = b + ((love.math.random() - 0.5) * n_chro + mono)

		return r, g, b, 1
	end)
	colour_map = love.graphics.newImage(id)
	vert_mesh_shader:send("colour_map", colour_map)
end
new_colours()

function draw_to_buffer()
	--
	lg.setShader(vert_mesh_shader)
	vert_mesh_shader:send("position_buffer", position_buffer)

	lg.setBlendMode("add", "premultiplied")
	
	--draw colour into screen buffer
	vert_mesh_shader:send("draw_colour", true)
	lg.setCanvas(screen_buffer)
	lg.draw(vert_mesh)

	--draw count
	vert_mesh_shader:send("draw_colour", false)
	lg.setCanvas(count_buffer)
	lg.draw(vert_mesh)

	lg.setCanvas()
	lg.setShader()
	lg.setBlendMode("alpha")
end

local total_iters = 0
local i = 0
local iters = 500
local trace_iters = 15
function _update()
	if i == 0 then
		random_positions()
		update_positions()
	else
		update_positions()
	end
	if i > trace_iters then
		draw_to_buffer()
	end
	i = i + 1
	if i > iters then
		i = 0
	end
	total_iters = total_iters + 1
end

function love.load()
	init_view()
end

function love.update(dt)
	update_cam(dt)

	local now = love.timer.getTime()
	while love.timer.getTime() - now < (10 / 1000) do
		_update()
	end
end

function love.draw()
	if view.dirty then
		lg.setCanvas(screen_buffer)
		lg.clear(0,0,0,0)
		lg.setCanvas(count_buffer)
		lg.clear(0,0,0,0)
		lg.setCanvas()
		view.dirty = false
		exposure = 1.0
		total_iters = 0
		i = 0
		for _ = 1, trace_iters * 2 do
			_update()
		end
	end
	
	local max_buffer = calc_max(count_buffer)
	screen_shader:send("maximum", max_buffer)
	screen_shader:send("count", count_buffer)

	lg.setShader(screen_shader)
	lg.draw(screen_buffer)
	lg.setShader()

	--debug
	if love.keyboard.isDown("`") then
		for i, v in ipairs {
			string.format("total iters: %5.3e - %d x %d (%d x %d)", total_iters * vert_count, total_iters, vert_count, vert_res, vert_res)
		} do
			lg.print(v, 10, 10 + (i - 1) * 20)
		end
		
		lg.draw(max_buffer, 10, 80)

		lg.draw(position_buffer, 10, 100)

		lg.draw(
			colour_map,
			0,
			lg.getHeight() - 2,
			0,
			lg.getWidth() / colour_map:getWidth(),
			2
		)
	end

end

function love.keypressed(k)
	local ctrl = love.keyboard.isDown("lctrl")
	if k == "q" and ctrl or k == "escape" then
		love.event.quit()
	end
	if k == "r" and ctrl then
		love.event.quit("restart")
	end

	if k == "m" then
		new_colours()
		view.dirty = true
	end

	if k == "n" then
		random_transforms()
		init_view()
	end
end
