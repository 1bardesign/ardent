local ffi = require("ffi")
require("lib.core")

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
	{
		{"VertexUV", "float", 2},
		{"VertexPosition", "float", 1}, --required to keep love happy
	},
	verts,
	"points",
	"static"
)
local vert_mesh_shader = lg.newShader([[
uniform Image position_buffer;
uniform Image colour_map;
uniform bool draw_colour;

uniform vec2 res;

varying float v_col;

#ifdef VERTEX
attribute vec2 VertexUV;
vec4 position(mat4 proj, vec4 _pos) {
	vec4 t = Texel(position_buffer, VertexUV);
	vec2 pos = t.rg;

	v_col = t.b;

	vec4 r = vec4(pos, 1.0, 1.0);

	return proj * r;
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

function init_view()
	view = {
		dirty = true,
		offset = vec2:zero(),
		scale = 200,
		angle = 0,
	}
end

function update_cam(dt)

	if not _g_old_view then
		_g_old_view = {}
		--copy the vectors
		for k,v in pairs(view) do
			if k ~= "dirty" then
				if type(v) == "table" then
					v = v:copy()
				end
				_g_old_view[k] = v
			end
		end
	end

	--translate
	local pan_speed = 500 / view.scale * dt
	local right = vec2:xy(1, 0):rotatei(-view.angle):smuli(pan_speed)
	local down = vec2:xy(0, 1):rotatei(-view.angle):smuli(pan_speed)
	if love.keyboard.isDown("w") then view.offset:vaddi(down:inverse()) end
	if love.keyboard.isDown("s") then view.offset:vaddi(down) end
	if love.keyboard.isDown("a") then view.offset:vaddi(right:inverse()) end
	if love.keyboard.isDown("d") then view.offset:vaddi(right) end

	--zoom
	local zoom_factor = 1 + dt
	if love.keyboard.isDown("z") then view.scale = view.scale * zoom_factor end
	if love.keyboard.isDown("x") then view.scale = view.scale / zoom_factor end

	--rotate
	local turn_speed = math.pi * 0.5 * dt
	if love.keyboard.isDown("q") then view.angle = view.angle + turn_speed end
	if love.keyboard.isDown("e") then view.angle = view.angle - turn_speed end

	--check for changes automagically
	local changed = false
	for k,v in pairs(_g_old_view) do
		local vk = view[k]
		if type(v) == "table" then
			if v:nequals(vk) then
				changed = true
			end
		elseif v ~= vk then
			changed = true
		end
		if changed then
			break
		end
	end
	if changed then
		view.dirty = true
		for k,v in pairs(_g_old_view) do
			local vk = view[k]
			if type(v) == "table" then
				v:vset(vk)
			else
				_g_old_view[k] = vk
			end
		end
	end
end

local transform_functions = {
	{"linear", ([[
		t = i;
	]]),},
	--multi component system
	{"sin", ([[
		t = sin(i);
	]]),},
	-- --single component wave
	{"sinx", ([[
		t = vec2(i.x, sin(i.x));
	]]),},
	--sqrt all components
	{"sqrt", ([[
		t = sign(i) * sqrt(abs(i));
	]]),},
	--classical flame transforms
	{"spherical", ([[
		t = i / (r2 + 0.01);
	]]),},
	{"swirl", ([[
		float c = cos(r);
		float s = sin(r);
		t = vec2(
			c * i.x - s * i.y,
			s * i.x + c * i.y
		);
	]]),},
	{"horseshoe", ([[
		t = vec2(
			(i.x - i.y) * (i.x + i.y),
			2.0 * i.x * i.y
		) / (r + 0.01);
	]]),},
	{"polar", ([[
		t = vec2(
			theta / pi,
			r - 1.0
		);
	]])},
	{"handkerchief", ([[
		t = r * vec2(
			sin(theta + r),
			cos(theta - r)
		);
	]])}

	--todo: distortion fns
	--noise
	--blur
}

--autogen cdef
ffi.cdef([[
struct flame_transform {
	//transforms
	float pre[6];
	float post[6];

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
		"\tvec2 t = vec2(0.0);",
		v[2], "\n",
		"\to = ", var_name, " * t + o;\n",
		"}",
	})
end
functions_snippet = table.concat(functions_snippet, "\n")

local position_buffer = lg.newCanvas(vert_res, vert_res, {format="rgba32f"})
local position_buffer_previous = lg.newCanvas(vert_res, vert_res, {format="rgba32f"})
local position_shader = lg.newShader([[
#pragma language glsl3

const int MAX_TRANSFORMS = 8;
const int TRANSFORM_SIZE = 6 + 6 + ]]..#transform_functions..[[;
uniform int transform_count;
uniform float transforms[MAX_TRANSFORMS * TRANSFORM_SIZE];

//random function sampling
uniform float transform_salt_offset;
uniform float transform_salt_scale;

//todo uniform
const float world_bounds = 100.0;

//
const float pi = 3.14159265359;

mat3 read_mat3(int idx) {
	return mat3(
		transforms[idx+0], transforms[idx+3], 0.0,
		transforms[idx+1], transforms[idx+4], 0.0,
		transforms[idx+2], transforms[idx+5], 1.0
	);
}

//from https://stackoverflow.com/questions/26070410/robust-atany-x-on-glsl-for-converting-xy-coordinate-to-angle
//more numerically stable
float atan2(float y, float x) {
	if (y == 0.0 && x == 0.0) {
		return 0.0;
	}
    float s = float(abs(x) > abs(y));
    return mix(pi / 2.0 - atan(x, y), atan(y, x), s);
}

vec2 tr(vec2 v_i, int tr) {
	int idx = tr * TRANSFORM_SIZE;

	///////////////////////////////////
	//collect transforms
	mat3 a = read_mat3(idx); idx += 6;
	mat3 b = read_mat3(idx); idx += 6;

	///////////////////////////////////
	//pre-transform
	vec2 i = (a * vec3(v_i, 1.0)).xy;

	///////////////////////////////////
	//apply activation

	//common parameters
	//distance from origin
	float r = length(i);
	float r2 = r * r;

	//angles
	float theta = atan2(i.y, i.x);
	float phi = atan2(i.x, i.y);

	//accumulate output in here
	vec2 o = vec2(0.0);

	]]..functions_snippet..[[

	///////////////////////////////////
	//post-transform
	o = (b * vec3(o, 1.0)).xy;

	//diverged? stay put
	//todo: re-init random?
	if(
		isnan(o.x) || isnan(o.y) ||
		isinf(o.x) || isinf(o.y)
	) {
		o = v_i;
	}

	o = sign(o) * min(abs(o), vec2(world_bounds));

	return o;
}

int get_transform(vec2 screenpos) {
	float px = floor(screenpos.x * transform_salt_scale) + transform_salt_offset;
	return (int(px) % transform_count);
}

#ifdef PIXEL
vec4 effect(vec4 color, Image tex, vec2 uv, vec2 screenpos) {
	float transform = get_transform(screenpos);
	//sample the point
	vec4 t = Texel(tex, uv);
	//apply
	return vec4(
		//new position
		tr(t.xy, int(transform)),
		//equidistant colours
		mix(t.z, float(transform) / (float(transform_count) - 1.0), 0.5),
		//constant alpha
		1.0
	);
}
#endif
]])

local transforms = {}

function upload_transforms()
	local count = #transforms
	local tr_data = love.data.newByteData(ffi.sizeof("struct flame_transform") * count)
	local tr_ser = ffi.cast("struct flame_transform*", tr_data:getFFIPointer())
	for i, transform in ipairs(transforms) do
		local ser = tr_ser[i - 1]

		--copy transforms
		for ti, v in ipairs(transform.pre) do
			ser.pre[ti - 1] = v
		end

		for ti, v in ipairs(transform.post) do
			ser.post[ti - 1] = v
		end

		for ti, v in ipairs(transform.function_weights) do
			ser.function_weights[ti - 1] = v
		end
	end
	position_shader:send("transforms", tr_data)
	position_shader:send("transform_count", count)
end

function chance(perc)
	return love.math.random() < perc
end

function random_transforms()
	--
	transforms = {}

	local function _random_rule()
		return love.math.random(1, #transform_functions)
	end
	--generate a certain number of transforms
	local num_transforms = love.math.random(3, 6)

	--amount of rules allowed per-transform
	local rules_max = 1 + math.floor(math.abs(love.math.randomNormal(0, 0.5)))

	--primary transform for the set
	local primary_rule = _random_rule()

	--minor rules
	local minor_rules_count = math.max(0, love.math.random(-2, 2))
	local minor_rules = {}
	while #minor_rules < minor_rules_count do
		table.insert(minor_rules, _random_rule())
	end
	local minor_rule_chance = minor_rules_count > 0 and love.math.random() * (1 / num_transforms) or 0

	--totally random rules
	local random_rules_allowed = chance(0.1)
	local random_rules_chance = random_rules_allowed and love.math.random() * (2 / num_transforms) or 0

	print("new config:")
	for i,v in ipairs {
		{"num_transforms", num_transforms},
		{"rules_max", rules_max},
		{"primary_rule", primary_rule},
		{"minor_rules_count", minor_rules_count, "chance", minor_rule_chance},
		{"random_rules_allowed", tostring(random_rules_allowed), "chance", random_rules_chance},
	} do
		print("", table.concat(v, "\t"))
	end

	function compute_transform(t)
		local c = math.cos(t.angle)
		local s = math.sin(t.angle)

		local px, py = t.offset_before:unpack()
		local ox, oy = t.offset_after:unpack()
		local sx, sy = t.scale:unpack()

		local rx = c * ox * sx - s * oy * sy + px
		local ry = s * ox * sx + c * oy * sy + py

		return {
			c * sx, -s * sy, rx,
			s * sx,  c * sy, ry,
		}
	end

	local function identity_transform()
		return {
			offset_before = vec2:xy(0, 0),
			offset_after = vec2:xy(0, 0),
			scale = vec2:xy(1, 1),
			angle = 0,
		}
	end

	local function mutate_transform(t, amount)
		--modify scale
		local function _r()
			return love.math.randomNormal(amount, 0)
		end
		local scale_modify_amount = 1 / 20
		if chance(0.85) then
			--uniform scale
			t.scale:smuli(
				1 - (amount * 0.5) + _r()
			)
		else
			--bilateral scale
			t.scale:smuli(
				1 - (amount * 0.5) + _r(),
				1 - (amount * 0.5) + _r()
			)
		end

		--modify angle
		t.angle = t.angle + _r() * 2 * math.pi

		t.offset_before:saddi(_r(), _r())

		if chance(0.5) then
			--dependent on offset before
			t.offset_after:vsubi(t.offset_before)
		end
		t.offset_after:saddi(_r(), _r())
		return t
	end

	local function rt()
		return mutate_transform(identity_transform(), math.lerp(0.1, 0.25, love.math.random()))
	end

	for i = 1, num_transforms do
		--randomly select a few rules with random weighting for this transform
		local rules_count  = love.math.random(1, rules_max)
		local function_weights = {}
		for i,v in ipairs(transform_functions_ordered) do function_weights[i] = 0 end
		for i = 0, rules_count do
			local f
			if chance(random_rules_chance) then
				f = _random_rule()
			elseif chance(minor_rules_count) then
				f = table.pick_random(minor_rules)
			else
				f = primary_rule
			end
			function_weights[f] = math.lerp(1, 5, love.math.random())
		end
		--normalise
		local t = 0
		for i,v in ipairs(function_weights) do
			v = math.abs(v)
			t = t + v
		end
		for i,v in ipairs(function_weights) do
			v = v / t
			function_weights[i] = v
		end

		local pre = rt()
		local post = rt()

		--output
		if true then
			print("\trule", i)
			for i,v in ipairs(function_weights) do
				if v ~= 0 then
					print("\t", transform_functions_ordered[i], v)
				end
			end
			for i,v in ipairs({
				{"pre", pre},
				{"post", post},
			}) do
				local name, t = unpack(v)
				print("\t", name, string.format(
					"scale: (% 5f % 5f) angle: % 5f offsets: (% 5f % 5f) (% 5f % 5f)",
					t.scale.x, t.scale.y,
					t.angle,
					t.offset_before.x, t.offset_before.y,
					t.offset_after.x, t.offset_after.y
				))
			end
		end

		local t = {
			pre = pre,
			post = post,
			function_weights = function_weights,
		}
		table.insert(transforms, t)
	end

	for i,v in ipairs(transforms) do
		v.pre = compute_transform(v.pre)
		v.post = compute_transform(v.post)
	end

	upload_transforms()
end

random_transforms()

function random_positions()
	position_buffer:setFilter("nearest", "nearest")
	position_buffer_previous:setFilter("nearest", "nearest")

	local id = love.image.newImageData(vert_res, vert_res, "rgba32f")
	id:mapPixel(function(x, y, r, g, b, a)
		x = (love.math.random() - 0.5) * 2
		y = (love.math.random() - 0.5) * 2
		return
			--pos
			x, y,
			--col
			love.math.random(),
			--unused
			1
	end)

	local im = love.graphics.newImage(id)
	lg.push("all")
	lg.setCanvas(position_buffer)
	lg.setBlendMode("replace", "premultiplied")
	lg.origin()
	lg.draw(im)
	lg.pop()

	upload_transforms()
end

function update_positions()
	--swap buffers
	position_buffer, position_buffer_previous = position_buffer_previous, position_buffer
	lg.setCanvas(position_buffer)
	lg.setShader(position_shader)
	lg.setBlendMode("replace", "premultiplied")
	--update random transform offset per-row
	do
		position_shader:send("transform_salt_offset", love.math.random(0, #transforms - 1))
		position_shader:send("transform_salt_scale", math.lerp(0.1, 0.5, love.math.random()))
	end

	--rotated to get a bit more scrambling across elements
	lg.draw(
		position_buffer_previous,
		vert_res, 0,
		math.pi * 0.5
	)

	--reset
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
screen_buffer = lg.newCanvas(
	sw, sh,
	{format="rgba32f"}
)
--uuuurgh we cant add on alpha channel currently (love 11.3 - fixed in 12.x apparently)
count_buffer = lg.newCanvas(
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
	//replace alpha from count buffer
	t.a = Texel(count, uv).r;

	//mean
	vec3 rgb = t.rgb / t.a;
	//linear
	//float alpha = t.a / highest;
	//loglog
	float alpha = log(1.0 + t.a) / log(1.0 + highest);

	//power / gamma
	alpha = pow(alpha, 1.0 / 2.5);

	//scale
	alpha *= 0.5;

	return vec4(
		rgb * alpha,
		alpha
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

function apply_camera()
	lg.origin()
	local sx, sy = screen_buffer:getDimensions()
	lg.translate(sx * 0.5, sy * 0.5)

	lg.scale(view.scale, view.scale)
	lg.rotate(view.angle)
	lg.translate(-view.offset.x, -view.offset.y)
end

function draw_to_buffer()
	--setup camera
	lg.push()
	apply_camera()

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

	lg.pop()
end

local total_iters = 0
local i = 0
local iters = 1000
local trace_iters = 20
function _update()
	if i == 0 then
		random_positions()
	end
	update_positions()

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
	local time_spin = 10 / 1000
	if love.keyboard.isDown("space") then
		time_spin = 1
	end
	while love.timer.getTime() - now < time_spin do
		_update()
	end
end

function love.draw()
	if view.dirty then
		view.dirty = false
		lg.setCanvas(screen_buffer)
		lg.clear(0,0,0,0)
		lg.setCanvas(count_buffer)
		lg.clear(0,0,0,0)
		lg.setCanvas()
		total_iters = 0
		i = 0
		for _ = 1, trace_iters * 2 do
			_update()
		end
	end

	local max_buffer = calc_max(count_buffer)
	screen_shader:send("maximum", max_buffer)
	screen_shader:send("count", count_buffer)

	lg.setBlendMode("alpha", "alphamultiply")

	--draw bg zero colour
	local sx, sy = screen_buffer:getDimensions()
	local cw, ch = colour_map:getDimensions()
	lg.setColor(1,1,1,0.5)
	lg.draw(
		colour_map, lg.newQuad(cw,0, 1,1, cw, ch),
		0, 0,
		0,
		sx, sy
	)
	lg.setColor(1,1,1,1)

	--draw screen
	lg.setBlendMode("alpha", "premultiplied")
	lg.setShader(screen_shader)
	lg.draw(screen_buffer)
	lg.setShader()
	lg.setBlendMode("alpha", "alphamultiply")

	--debug
	if love.keyboard.isDown("`") then
		for i, v in ipairs {
			string.format("total iters: %5.3e - %d x %d (%d x %d)", total_iters * vert_count, total_iters, vert_count, vert_res, vert_res)
		} do
			lg.print(v, 10, 10 + (i - 1) * 20)
		end

		lg.draw(max_buffer, 10, 80)

		if not position_buffer_shader then
			position_buffer_shader = lg.newShader(nil, [[
				vec4 effect(vec4 col, Image tex, vec2 uv, vec2 screenpos) {
					return vec4(vec3(0.5) + Texel(tex, uv).rgb * 0.25, 1.0);
				}
			]])
		end
		lg.setShader(position_buffer_shader)
		lg.draw(position_buffer, 10, 100)
		lg.setShader()

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

	if k == "m" or k == "n" then
		new_colours()
		view.dirty = true
	end

	if k == "n" or k == "j" then
		random_transforms()
		init_view()
	end
end
