--matrix stuff

--internal matrix multiplication element
local function _mm(a, b, i, j, aw, ah, bw, bh)
	local r = 0
	i = i - 1
	j = j - 1
	for c = 0, aw - 1 do
		local ai = (c * ah + i) + 1
		local bi = (j * bw + c) + 1
		r = r + a[ai] * b[bi]
	end
	return r
end

--generic matmul
--returns the matrix, and its dimensions
function matmul(a, aw, ah, b, bw, bh)
	if aw ~= bh then
		error("matrix multiplication with mismatched matrices")
	end
	local r = {}
	for i = 1, ah do
		for j = 1, bw do
			local v = _mm(a, b, i, j, aw, ah, bw, bh)
			table.insert(r, v)
		end
	end
	return r, bw, ah
end

function matmul4x4(a, b)
	return matmul(a, 4, 4, b, 4, 4)
end

function matmul_v4(v, m)
	return matmul(v, 4, 1, m, 4, 4)
end

function transform_point(v, m)
	local res = matmul_v4({v.x, v.y, v.z, 1.0}, m)
	return vec3:xyz(res[1], res[2], res[3])
end

function transform_direction(v, m)
	local res = matmul_v4({v.x, v.y, v.z, 0.0}, m)
	return vec3:xyz(res[1], res[2], res[3])
end

--3d view stuff

function calculate_perspective(screen, fov, near, far, into)
	--everything we need
	local aspect = screen[1] / screen[2];
	local f = 1.0 / math.tan(fov / 2.0);
	local xpr = f / aspect;
	local ypr = f;
	local fmn = (far - near);
	local zpr = (far + near) / fmn;
	local zhpr = (2.0 * far * near) / fmn;
	
	--garbage generation
	into = into or {}

	--zero out
	for i = 1, 16 do
		into[i] = 0
	end
	--fill required elements
	--(row major, opposite of in glsl)
	into[1] = xpr
	into[6] = ypr
	into[11] = zpr
	into[12] = zhpr
	into[15] = 1

	return into
end

function calculate_camera_core(right, up, forward, from)
	--todo: ensure ortho
	local d_r = -from:dot(right)
	local d_u = -from:dot(up)
	local d_f = -from:dot(forward)
	return {
		right.x,   right.y,   right.z,   d_r,
		up.x,      up.y,      up.z,      d_u,
		forward.x, forward.y, forward.z, d_f,
		0.0,       0.0,       0.0,       1.0,
	}
end

function calculate_camera(from, to, natural_up)
	local forward = from:vsub(to):normalisei()
	local right = natural_up:cross(forward):normalisei()
	local up = forward:cross(right):normalisei()
	return calculate_camera_core(right, up, forward, from)
end

function rotation_axis_angle(v, a)
	--via iq, copied with appropriate shame :)
    local si = math.sin(a)
    local co = math.cos(a)
    local ic = 1.0 - co

    return {
    	v.x*v.x*ic + co,       v.y*v.x*ic - si*v.z,    v.z*v.x*ic + si*v.y, 0,
        v.x*v.y*ic + si*v.z,   v.y*v.y*ic + co,        v.z*v.y*ic - si*v.x, 0,
        v.x*v.z*ic - si*v.y,   v.y*v.z*ic + si*v.x,    v.z*v.z*ic + co,     0,
        0,                     0,                      0,                   1,
    }
end
