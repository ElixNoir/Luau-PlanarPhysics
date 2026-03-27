--- Author: COMT

--[[

	PlanarPhysics - A module for handling planar physics.
	Mod: RigidBody
	Version: 0.0.1

]]


--!native
--!optimize 2
--!strict


-- -- State

local Manifold = require("./Manifold")
type Manifold = Manifold.Manifold

local Module = {}
Module.__index = Module

export type RigidBody = {
	Angle: number,
	AngularVelocity: number,
	Friction_Kinetic: number, Friction_Static: number,
	Inertia: number, Inertia_Inverse: number,
	Mass: number, Mass_Inverse: number,
	Position_X: number, Position_Y: number,
	Shape: {number},
	Velocity_X: number, Velocity_Y: number,
	
	update_AABB: (self: RigidBody) -> (),

	-- Internal
	_AABB_Minimum_X: number, _AABB_Minimum_Y: number,
	_AABB_Maximum_X: number, _AABB_Maximum_Y: number,
	_Island_Identity: number,
	_Sleeping: boolean,
	_Sleep_Timer: number,
	_Index: number -- position in world bodelta_Y list
}


-- -- Helpers

local math_abs = math.abs
local math_clamp = math.clamp
local math_cos = math.cos
local math_Huge = math.huge
local math_max = math.max
local math_min = math.min
local math_sin = math.sin
local math_sqrt = math.sqrt
local math_Pi = math.pi

local table_clear = table.clear
local table_create = table.create
local table_sort = table.sort


-- -- Functions

-- Applocal_Ya force (accumulates into velocity for one step).
function Module.apply_Force(self: RigidBody, fx: number, fy: number, deltaTime: number)
	self.Velocity_X += fx * self.Mass_Inverse * deltaTime
	self.Velocity_Y += fy * self.Mass_Inverse * deltaTime
	self._Sleeping = false
	self._Sleep_Timer = 0
end

-- Applocal_Yan impulse (world space) at a world-space point.
function Module.apply_Impulse(self: RigidBody, jx: number, jy: number, px: number, py: number)
	local rx, ry = px - self.Position_X, py - self.Position_Y
	self.AngularVelocity += (rx * jy - ry * jx) * self.Inertia_Inverse
	self.Velocity_X += jx * self.Mass_Inverse
	self.Velocity_Y += jy * self.Mass_Inverse
	self._Sleeping = false
	self._Sleep_Timer = 0
end


-- -- -- Circle

function Module.new_Circle(radius: number): RigidBody
	return setmetatable({
		Angle = 0,
		AngularVelocity = 0,
		Friction_Static = .5, Friction_Kinetic = .3,
		Inertia = 1, Inertia_Inverse = 1,
		Mass = 1, Mass_Inverse = 1,
		Position_X = 0, Position_Y = 0,
		Shape = {radius},
		Velocity_X = 0, Velocity_Y = 0,
		
		update_AABB = Module.update_Circle_AABB,

		_AABB_Minimum_X = 0, _AABB_Minimum_Y = 0,
		_AABB_Maximum_X = 0, _AABB_Maximum_Y = 0,
		_Island_Identity = 0,
		_Sleeping = false,
		_Sleep_Timer = 0,
		_Index = 0
	}, Module)
end

function Module.collide_Circle_Circle(circleA: RigidBody, circleB: RigidBody, manifold: Manifold): boolean
	local delta_X, delta_Y = circleB.Position_X - circleA.Position_X, circleB.Position_Y - circleA.Position_Y
	local distanceSquared = delta_X * delta_X + delta_Y * delta_Y
	
	local a_Radius = circleA.Shape[1]
	local combinedRadius = a_Radius + circleB.Shape[1]
	
	if distanceSquared >= combinedRadius * combinedRadius then
		return false
	end

	local manifold_Contact_Normal_X: number, manifold_Contact_Normal_Y: number

	local distance = math_sqrt(distanceSquared)
	if distance < 1e-9 then
		manifold_Contact_Normal_X, manifold_Contact_Normal_Y = 1, 0
		manifold.Penetration = combinedRadius
	else
		local inverseDistance = 1 / distance
		manifold_Contact_Normal_X, manifold_Contact_Normal_Y = delta_X * inverseDistance, delta_Y * inverseDistance
		manifold.Penetration = combinedRadius - distance
	end

	-- Single contact point at surface of A
	manifold.Contact1_X = circleA.Position_X + manifold_Contact_Normal_X * a_Radius
	manifold.Contact1_Y = circleA.Position_Y + manifold_Contact_Normal_Y * a_Radius
	manifold.Contact_Count = 1

	manifold.Contact_Normal_X, manifold.Contact_Normal_Y = manifold_Contact_Normal_X, manifold_Contact_Normal_Y

	return true
end

function Module.collide_Circle_OBB(circle: RigidBody, box: RigidBody, manifold: Manifold, circle_is_a: boolean): boolean
	-- Express circle centre in box local space
	local box_Angle = box.Angle
	local delta_X, delta_Y = circle.Position_X - box.Position_X, circle.Position_Y - box.Position_Y
	local cosine, sine = math_cos(box_Angle), math_sin(-box_Angle)
	local local_X, local_Y= cosine * delta_X - sine * delta_Y, sine * delta_X + cosine * delta_Y

	-- Find half-extents from shape vertices (flat [x1,y1,…])
	local shape = box.Shape
	local hx, hy = 0, 0
	for index = 1, #shape, 2 do
		local ax, ay = math_abs(shape[index]), math_abs(shape[index + 1])
		if ax > hx then hx = ax end
		if ay > hy then hy = ay end
	end

	-- Clamp circle centre to box
	local cx_clamped = math_clamp(local_X, -hx, hx)
	local cy_clamped = math_clamp(local_Y, -hy, hy)
	local inside = (local_X == cx_clamped) and (local_Y== cy_clamped)

	local fx = local_X - cx_clamped
	local fy = local_Y- cy_clamped
	local distanceSquared = fx * fx + fy * fy
	local radius = circle.Shape[1]

	if not inside and distanceSquared >= radius * radius then
		return false
	end

	local nx_local, ny_local: number, number
	local penetration: number

	if inside then
		-- Push out along nearest face
		local overlap_X, overlap_Y = hx - math_abs(local_X), hy - math_abs(local_Y)
		if overlap_X < overlap_Y then
			penetration = radius + overlap_X
			nx_local = if local_X > 0 then 1 else -1
			ny_local = 0
		else
			penetration = radius + overlap_Y
			nx_local = 0
			ny_local = if local_Y> 0 then 1 else -1
		end
	else
		local distance = math_sqrt(distanceSquared)
		penetration = radius - distance
		local inverseDistance = 1 / distance
		nx_local = fx * inverseDistance
		ny_local = fy * inverseDistance
	end

	-- Rotate normal back to world space

	local nx_world, ny_world = cosine * nx_local + sine * ny_local, -sine * nx_local + cosine * ny_local

	local manifold_Contact_Normal_X: number, manifold_Contact_Normal_Y: number
	if circle_is_a then
		manifold_Contact_Normal_X, manifold_Contact_Normal_Y = -nx_world, -ny_world
	else
		manifold_Contact_Normal_X, manifold_Contact_Normal_Y = nx_world, ny_world
	end

	manifold.Contact_Normal_X, manifold.Contact_Normal_Y = manifold_Contact_Normal_X, manifold_Contact_Normal_Y

	manifold.Penetration = penetration
	manifold.Contact1_X = circle.Position_X + manifold_Contact_Normal_X * radius
	manifold.Contact1_Y = circle.Position_Y + manifold_Contact_Normal_Y * radius
	manifold.Contact_Count = 1

	return true
end

function Module.update_Circle_AABB(self: RigidBody)
	local position_X, position_Y = self.Position_X, self.Position_Y
	local radius = self.Shape[1]
	self._AABB_Maximum_X, self._AABB_Maximum_Y = position_X + radius, position_Y + radius
end


-- -- -- Polygon

function Module.new_Polygon(shape: {number}): RigidBody
	return setmetatable({
		Angle = 0,
		AngularVelocity = 0,
		Friction_Static = .5, Friction_Kinetic = .3,
		Inertia = 1, Inertia_Inverse = 1,
		Mass = 1, Mass_Inverse = 1,
		Position_X = 0, Position_Y = 0,
		Shape = shape,
		Velocity_X = 0, Velocity_Y = 0,

		update_AABB = Module.update_Polygon_AABB,

		_AABB_Minimum_X = 0, _AABB_Minimum_Y = 0,
		_AABB_Maximum_X = 0, _AABB_Maximum_Y = 0,
		_Island_Identity = 0,
		_Sleeping = false,
		_Sleep_Timer = 0,
		_Index = 0
	}, Module)
end

-- Project polygon onto axis, return [min,max]
local function project_poly(verts: { number }, ax: number, ay: number): (number, number)
	local min_p, max_p = math_Huge, -math_Huge

	for index = 1, #verts, 2 do
		local x, y = verts[index], verts[index + 1]
		local dot = x * ax + y * ay
		if dot < min_p then min_p = dot end
		if dot > max_p then max_p = dot end
	end

	return min_p, max_p
end

-- Find best (least-penetrating) axis from one polygon against another.
-- Returns: overlap, contact_Normal_X, contact_Normal_Y, edge_vertex_index
local function best_axis(a_verts: { number }, b_verts: { number }): (number, number, number, number)
	local length = #a_verts
	local best_sep = -math_Huge
	local best_nx = 0
	local best_ny = 0
	local best_edge = 1
	for index = 1, length, 2 do
		local j = (index % length) + 1
		-- Compute edge normal
		local ex = a_verts[index + 2] - a_verts[index]
		-- Safe pair indexing
		local cur_x = a_verts[index]
		local cur_y = a_verts[index + 1]
		local ni = if index + 2 <= length then index + 2 else 1
		local next_x = a_verts[ni]
		local next_y = a_verts[ni + 1]
		local edge_x = next_x - cur_x
		local edge_y = next_y - cur_y
		-- Outward normal (right-hand)
		local len = math_sqrt(edge_x * edge_x + edge_y * edge_y)
		if len < 1e-10 then
			continue
		end
		local contact_Normal_X = edge_y / len
		local contact_Normal_Y = -edge_x / len

		local a_min, a_max = project_poly(a_verts, contact_Normal_X, contact_Normal_Y)
		local b_min, b_max = project_poly(b_verts, contact_Normal_X, contact_Normal_Y)
		local sep = math_max(a_min - b_max, b_min - a_max)

		if sep > 0 then
			return sep, contact_Normal_X, contact_Normal_Y, index
		end
		if sep > best_sep then
			best_sep = sep
			best_nx = contact_Normal_X
			best_ny = contact_Normal_Y
			best_edge = index
		end
	end
	return best_sep, best_nx, best_ny, best_edge
end

-- Clip an edge (incident) against a side plane defined by (ref_v, ref_n).
local function clip_segment(
	p1x: number, p1y: number, p2x: number, p2y: number,
	contact_Normal_X: number, contact_Normal_Y: number, offset: number
): (number, number, number, number, number)
	-- Returns clipped p1x,p1y, p2x,p2y, count
	local d1, d2 = p1x * contact_Normal_X + p1y * contact_Normal_Y - offset, p2x * contact_Normal_X + p2y * contact_Normal_Y - offset
	if d1 >= 0 and d2 >= 0 then
		return p1x, p1y, p2x, p2y, 2
	end
	if d1 < 0 and d2 < 0 then
		return p1x, p1y, p2x, p2y, 0
	end
	local t = d1 / (d1 - d2)
	local cx = p1x + t * (p2x - p1x)
	local cy = p1y + t * (p2y - p1y)
	if d1 < 0 then
		return cx, cy, p2x, p2y, 2
	end
	return p1x, p1y, cx, cy, 2
end

-- Returns world-space vertices of a polygon bodelta_Y
local function world_verts(self: RigidBody): { number }
	local angle = self.Angle
	local position_X, position_Y = self.Position_X, self.Position_Y
	local shape = self.Shape

	local cosine, sine = math_cos(angle), math_sin(angle)

	local out: { number } = {}
	for index = 1, #shape, 2 do
		local local_X, local_Y= shape[index], shape[index + 1]
		out[#out + 1], out[#out + 2] = cosine * local_X - sine * local_Y+ position_X, sine * local_X + cosine * local_Y+ position_Y
	end

	return out
end


-- -- Methods

function Module.collide_Polygon(a: RigidBody, b: RigidBody, manifold: Manifold): boolean
	local va = world_verts(a)
	local vb = world_verts(b)

	local sep_a, nx_a, ny_a = best_axis(va, vb)
	if sep_a > 0 then
		return false
	end
	
	local sep_b, nx_b, ny_b = best_axis(vb, va)
	if sep_b > 0 then
		return false
	end

	local contact_Normal_X, contact_Normal_Y: number, number
	if sep_a > sep_b then
		contact_Normal_X, contact_Normal_Y = nx_a, ny_a
	else
		contact_Normal_X, contact_Normal_Y = -nx_b, -ny_b
	end

	-- Ensure normal points A→B
	local delta_X = b.Position_X - a.Position_X
	local delta_Y = b.Position_Y - a.Position_Y
	if delta_X * contact_Normal_X + delta_Y * contact_Normal_Y < 0 then
		contact_Normal_X = -contact_Normal_X
		contact_Normal_Y = -contact_Normal_Y
	end

	-- Find incident edge on B (most anti-parallel to normal)
	local nb = #vb
	local best_dot = math_Huge
	local inc_i = 1
	local ii = 1
	while ii < nb do
		local ni2 = if ii + 2 <= nb then ii + 2 else 1
		local ex = vb[ni2] - vb[ii]
		local ey = vb[ni2 + 1] - vb[ii + 1]
		local len = math_sqrt(ex * ex + ey * ey)
		local q1, q2 = ex / len, ey / len
		local d = q1 * contact_Normal_X + q2 * contact_Normal_Y
		if d < best_dot then
			best_dot = d
			inc_i = ii
		end
		ii += 2
	end
	local ni2 = if inc_i + 2 <= nb then inc_i + 2 else 1
	local ip1x, ip1y = vb[inc_i], vb[inc_i + 1]
	local ip2x, ip2y = vb[ni2], vb[ni2 + 1]

	-- Side planes from reference face
	local tang_x = contact_Normal_Y; local tang_y = -contact_Normal_X -- tangent along reference edge
	local na = #va
	local ref_i = 1
	local best_ref = -math_Huge
	local ri = 1
	while ri < na do
		local q1, q2 = va[ri], va[ri + 1]
		local d = q1 * contact_Normal_X + q2 * contact_Normal_Y
		if d > best_ref then best_ref = d; ref_i = ri end
		ri += 2
	end
	local ref_x, ref_y = va[ref_i], va[ref_i + 1]
	local offset1 = ref_x * tang_x + ref_y * tang_y
	local offset2 = ref_x * tang_x + ref_y * tang_y -- opposite side

	local c1, count1: number
	ip1x, ip1y, ip2x, ip2y, count1 = clip_segment(ip1x, ip1y, ip2x, ip2y, tang_x, tang_y, offset1)
	if count1 < 2 then
		return false
	end
	ip1x, ip1y, ip2x, ip2y, count1 = clip_segment(ip1x, ip1y, ip2x, ip2y, -tang_x, -tang_y, offset2)
	if count1 < 2 then
		return false
	end

	-- Keep points below reference face
	local ref_dot = ref_x * contact_Normal_X + ref_y * contact_Normal_Y
	local pen1 = ref_dot - ip1x * contact_Normal_X - ip1y * contact_Normal_Y
	local pen2 = ref_dot - ip2x * contact_Normal_X - ip2y * contact_Normal_Y
	local count = 0

	manifold.Contact_Count = 0
	if pen1 >= 0 then
		manifold.Contact1_X, manifold.Contact1_Y = ip1x, ip1y
		count += 1
	end
	if pen2 >= 0 then
		if count == 0 then
			manifold.Contact1_X, manifold.Contact1_Y = ip2x, ip2y
		else
			manifold.Contact2_X, manifold.Contact2_Y = ip2x, ip2y
		end
		count += 1
	end
	if count == 0 then
		return false
	end

	manifold.Contact_Count = count
	manifold.Penetration = if pen1 > pen2 then pen1 else pen2
	manifold.Contact_Normal_X = contact_Normal_X
	manifold.Contact_Normal_Y = contact_Normal_Y

	return true
end

function Module.update_Polygon_AABB(self: RigidBody)
	local angle = self.Angle
	local shape = self.Shape
	local position_X, position_Y = self.Position_X, self.Position_Y

	local cosine, sine = math_cos(angle), math_sin(angle)
	local min_x, min_y = math_Huge, math_Huge
	local max_x, max_y = -math_Huge, -math_Huge

	for index = 1, #shape, 2 do
		local local_X, local_Y= shape[index], shape[index + 1]

		local wx, wy = cosine * local_X - sine * local_Y + position_X, sine * local_X + cosine * local_Y + position_Y
		if wx < min_x then min_x = wx end
		if wy < min_y then min_y = wy end
		if wx > max_x then max_x = wx end
		if wy > max_y then max_y = wy end
	end

	self._AABB_Maximum_X, self._AABB_Maximum_Y = max_x, max_y
end


-- -- End

return Module
