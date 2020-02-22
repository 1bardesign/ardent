local ffi = require("ffi")
lg = love.graphics

local vert_res = 50
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

varying float v_col;

#ifdef VERTEX
attribute vec2 VertexUV;
vec4 position(mat4 proj, vec4 _pos) {
	vec4 t = Texel(position_buffer, VertexUV);
	vec2 pos = t.xy;
	v_col = t.z;
	return proj * vec4(pos, 0.0, 1.0);
}
#endif
#ifdef PIXEL
vec4 effect(vec4 col, Image tex, vec2 uv, vec2 screenpos) {
	if (!draw_colour) {
		return vec4(1.0);
	}
	return Texel(colour_map, vec2(v_col, 0.0));
}
#endif
]])

ffi.cdef[[
struct flame_transform{
	//pre-transform
	float axx, axy, axc;
	float ayx, ayy, ayc;

	//post-transform
	float bxx, bxy, bxc;
	float byx, byy, byc;

	//function lerps
	float f_linear;
	float f_sin;
	float f_cos;
	float f_tan;
	float f_sqrt;
	float f_spherical;
};
]]

--todo: autogen cdef

local transform_functions = {
	{"linear", 		"i",},
	{"sin", 		"sin(i)",},
	{"cos", 		"cos(i)",},
	{"tan", 		"tan(i)",},
	{"sqrt", 		"sign(i) * sqrt(abs(i))",},
	{"spherical", 	"i / (r * r)",},
}

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
		"\to += ", var_name, " * ", v[2], ";\n",
		"}",
	})
end
functions_snippet = table.concat(functions_snippet, "\n")

local position_buffer = lg.newCanvas(vert_res, vert_res, {format="rgba32f"})
local position_buffer_previous = lg.newCanvas(vert_res, vert_res, {format="rgba32f"})
local position_shader = lg.newShader([[
const int MAX_TRANSFORMS = 8;
const int TRANSFORM_SIZE = 6 + 6 + 6;
uniform int transform_count;
uniform int transform_offset;
uniform int transform_stride;
uniform float transforms[MAX_TRANSFORMS * TRANSFORM_SIZE];

vec2 tr(vec2 i, int tr) {
	int idx = tr * TRANSFORM_SIZE;

	///////////////////////////////////
	//pre-transform
	float axx = transforms[idx++];
	float axy = transforms[idx++];
	float axc = transforms[idx++];

	float ayx = transforms[idx++];
	float ayy = transforms[idx++];
	float ayc = transforms[idx++];

	i = vec2(
		axx * i.x + axy * i.y + axc,
		ayx * i.x + ayy * i.y + ayc
	);

	///////////////////////////////////
	//post-transform
	float bxx = transforms[idx++];
	float bxy = transforms[idx++];
	float bxc = transforms[idx++];

	float byx = transforms[idx++];
	float byy = transforms[idx++];
	float byc = transforms[idx++];

	///////////////////////////////////
	//apply activation
	float r = length(i);

	vec2 o = vec2(0.0);
	
	]]..functions_snippet..[[

	//(apply post)
	o = vec2(
		bxx * o.x + bxy * o.y + bxc,
		byx * o.x + byy * o.y + byc
	);

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
		tr(t.xy, int(transform)),
		//equidistant colours
		mix(t.z, (float(transform) + 0.5) / float(transform_count), 0.5),
		//full alpha
		1.0
	);
}
#endif
]])

local transforms = {
	{
		pre = {
			offset = {0.3, 0.1},
			angle = 0.5,
			around = {1.0, 0.0},
			scale = {0.95, 0.95},
		},
		post = {
			offset = {0.0, 0.0},
			angle = 0.0,
			around = {0.0, 0.0},
			scale = {1.0, 1.0},
		},
		func = "sqrt",
	},
	{
		pre = {
			offset = {0.0, -1.0},
			angle = 0.3,
			around = {-1.0, 1.0},
			scale = {1.1, 1.1},
		},
		post = {
			offset = {0.0, 0.0},
			angle = 0.0,
			around = {0.0, 0.0},
			scale = {1.0, 1.0},
		},
		func = "spherical",
	},
	{
		pre = {
			offset = {1.0, 0.0},
			angle = -0.2,
			around = {0.0, 1.0},
			scale = {0.95, 0.95},
		},
		post = {
			offset = {0.0, 0.0},
			angle = 0.0,
			around = {0.0, 0.0},
			scale = {1.0, 1.0},
		},
		func = "spherical",
	},
	{
		pre = {
			offset = {1.3, 1.3},
			angle = -0.7,
			around = {1.0, 1.0},
			scale = {1.15, 1.15},
		},
		post = {
			offset = {0.0, 0.0},
			angle = 0.0,
			around = {0.0, 0.0},
			scale = {1.0, 1.0},
		},
		func = "spherical",
	},
}

function to_tr(t)
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

		--transforms
		t.axx, t.axy, t.axc, t.ayx, t.ayy, t.ayc = unpack(to_tr(v.pre))
		t.bxx, t.bxy, t.bxc, t.byx, t.byy, t.byc = unpack(to_tr(v.post))

		--normalised rule amounts
		local rule_amounts = {}
		local func = v.func
		for i, v in ipairs(transform_functions_ordered) do
			local f = 0
			if type(func) == "table" then
				if func[v] then
					--set table
					f = func[v]
				else
					--ordered table of either {"f", v} or "f"
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
		t.f_linear =    rule_amounts[1]
		t.f_sin =       rule_amounts[2]
		t.f_cos =       rule_amounts[3]
		t.f_tan =       rule_amounts[4]
		t.f_sqrt =      rule_amounts[5]
		t.f_spherical = rule_amounts[6]
	end
	position_shader:send("transforms", tr_d)
	position_shader:send("transform_count", count)
end

function random_transforms()
	--
	transforms = {}
	for i = 1, love.math.random(2, 8) do
		--randomly select 1-3 rules with random weighting
		local function_weights = {}
		for i,v in ipairs(transform_functions_ordered) do
			function_weights[i] = 0
		end
		for i = 1, love.math.random(1, 3) do
			function_weights[love.math.random(1, #function_weights)] = 0.25 + love.math.random()
		end
		--normalise
		local t = 0
		for i,v in ipairs(function_weights) do
			t = t + function_weights[i]
		end
		for i,v in ipairs(function_weights) do
			function_weights[i] = v / t
		end

		local function rt()
			return {
				offset = {
					love.math.randomNormal(),
					love.math.randomNormal(),
				},
				angle = (love.math.random() * 2 - 1) * math.pi,
				around = {
					love.math.randomNormal(),
					love.math.randomNormal(),
				},
				scale = {
					love.math.randomNormal(1),
					love.math.randomNormal(1),
				},
			}
		end

		local t = {
			pre = rt(),
			post = rt(),
			func = function_weights,
		}
		table.insert(transforms, t)
	end
	upload_transforms()
end

function random_positions()
	position_buffer:setFilter("nearest", "nearest")
	position_buffer_previous:setFilter("nearest", "nearest")

	lg.setCanvas(position_buffer)
	lg.clear(0,0,0,0)
	for i = 1, 8 do
		lg.setBlendMode(
			(i % 2) == 0 and "add"
			or "subtract"
		)
		for y = 1, vert_res do
			for x = 1, vert_res do
				lg.setColor(
					love.math.random(),
					love.math.random(),
					love.math.random(),
					1
				)
				lg.points(x - 1, y - 1)
			end
		end
	end
	lg.setColor(1,1,1,1)
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

	--rotated to ensure proper scrambling
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

local view = {
	dirty = true,
	offset = {0, 0},
	zoom = 100,
	rotation = 0,
}

function draw_to_buffer()
	lg.push()
	lg.origin()
	local bw, bh = screen_buffer:getDimensions()
	
	lg.translate(bw * 0.5, bh * 0.5)

	lg.translate(view.offset[1], view.offset[2])
	lg.rotate(view.rotation)
	lg.scale(view.zoom, view.zoom)

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
	lg.pop()
end

local total_iters = 0
local i = 0
local iters = 100
local trace_iters = 10
function _update()
	if i == 0 then
		random_positions()
		for _ = 1, trace_iters do
			update_positions()
		end
	else
		update_positions()
	end
	draw_to_buffer()
	i = i + 1
	if i > iters then
		i = 0
	end
	total_iters = total_iters + 1
end

function love.load()
end

function love.update(dt)
	local now = love.timer.getTime()
	while love.timer.getTime() - now < (10 / 1000) do
		_update()
	end

	local vo = view.offset
	local move_step = dt * 100
	local zoom_scale = 1 + 0.5 * dt
	local rotation_step = dt * math.pi * 0.1

	if love.keyboard.isDown("w") then view.dirty = true; vo[2] = vo[2] + move_step end
	if love.keyboard.isDown("s") then view.dirty = true; vo[2] = vo[2] - move_step end
	if love.keyboard.isDown("a") then view.dirty = true; vo[1] = vo[1] + move_step end
	if love.keyboard.isDown("d") then view.dirty = true; vo[1] = vo[1] - move_step end
	if love.keyboard.isDown("q") then view.dirty = true; view.zoom = view.zoom / zoom_scale end 
	if love.keyboard.isDown("e") then view.dirty = true; view.zoom = view.zoom * zoom_scale end 
	if love.keyboard.isDown("r") then view.dirty = true; view.rotation = view.rotation - rotation_step end
	if love.keyboard.isDown("f") then view.dirty = true; view.rotation = view.rotation + rotation_step end
	if love.keyboard.isDown("c") then view.dirty = true end
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
		draw_to_buffer()
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
			string.format("total iters: %d x %d (%d x %d)", total_iters, vert_count, vert_res, vert_res)
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
		view.dirty = true
	end
end
