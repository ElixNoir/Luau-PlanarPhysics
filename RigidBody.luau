--- Author: COMT

--[[

	PlanarPhysics - A module for handling planar physics.
	Mod: RigidBody
	Version: 0.0.2

]]


--!native
--!optimize 2
--!strict


-- -- State

local Manifold = require("./Manifold")
type Manifold = Manifold.Manifold

local Module = {}
Module.__index = Module

export type RigidBody = typeof(setmetatable({} :: {
	Angle: number,
	AngularVelocity: number,
	Friction_Kinetic: number, Friction_Static: number,
	Inertia: number, Inertia_Inverse: number,
	Mass: number, Mass_Inverse: number,
	Position_X: number, Position_Y: number,
	Vertices: {number},
	Velocity_X: number, Velocity_Y: number,
	
	update_AABB: (self: RigidBody) -> (),

	-- Internal
	_AABB_Minimum_X: number, _AABB_Minimum_Y: number,
	_AABB_Maximum_X: number, _AABB_Maximum_Y: number,
	_Island_Identity: number,
	_Sleeping: boolean,
	_Sleep_Timer: number,
	_Index: number
}, Module))


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

--[[
	Apply a force to the RigidBody.
]]
function Module.apply_Force(self: RigidBody, force_X: number, force_Y: number, deltaTime: number)
	self.Velocity_X += force_X * self.Mass_Inverse * deltaTime
	self.Velocity_Y += force_Y * self.Mass_Inverse * deltaTime
	self._Sleeping = false
	self._Sleep_Timer = 0
end

--[[
	Apply an impulse to the RigidBody at a given position.
]]
function Module.apply_Impulse(self: RigidBody, impulse_X: number, impulse_Y: number, position_X: number, position_Y: number)
	local rx, ry = position_X - self.Position_X, position_Y - self.Position_Y
	self.AngularVelocity += (rx * impulse_Y - ry * impulse_X) * self.Inertia_Inverse
	self.Velocity_X += impulse_X * self.Mass_Inverse
	self.Velocity_Y += impulse_Y * self.Mass_Inverse
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
		Vertices = {radius},
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

--[[
	Check collision between two circles.
	returns: collided
]]
function Module.collide_Circle_Circle(self: RigidBody, circle: RigidBody, manifold: Manifold): boolean
	local delta_X, delta_Y = circle.Position_X - self.Position_X, circle.Position_Y - self.Position_Y
	local distanceSquared = delta_X * delta_X + delta_Y * delta_Y
	
	local circleA_Radius = self.Vertices[1]
	local combinedRadius = circleA_Radius + circle.Vertices[1]
	
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

	manifold.Contact1_X = self.Position_X + manifold_Contact_Normal_X * circleA_Radius
	manifold.Contact1_Y = self.Position_Y + manifold_Contact_Normal_Y * circleA_Radius
	manifold.Contact_Count = 1
	manifold.Contact_Normal_X, manifold.Contact_Normal_Y = manifold_Contact_Normal_X, manifold_Contact_Normal_Y

	return true
end

--[[
	Check collision between a circle and a polygon.
	returns: collided
]]
function Module.collide_Circle_Polygon(self: RigidBody, polygon: RigidBody, manifold: Manifold, circle_is_a: boolean): boolean
	-- Express circle centre in box local space
	local box_Angle = polygon.Angle
	local delta_X, delta_Y = self.Position_X - polygon.Position_X, self.Position_Y - polygon.Position_Y
	local cosine, sine = math_cos(box_Angle), math_sin(-box_Angle)
	local local_X, local_Y= cosine * delta_X - sine * delta_Y, sine * delta_X + cosine * delta_Y

	-- Find half-extents from vertices (flat [x1,y1,…])
	local vertices = polygon.Vertices
	local hx, hy = 0, 0
	for index = 1, #vertices, 2 do
		local ax, ay = math_abs(vertices[index]), math_abs(vertices[index + 1])
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
	local radius = self.Vertices[1]

	if not inside and distanceSquared >= radius * radius then
		return false
	end

	local nx_local: number, ny_local: number
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
	manifold.Contact1_X = self.Position_X + manifold_Contact_Normal_X * radius
	manifold.Contact1_Y = self.Position_Y + manifold_Contact_Normal_Y * radius
	manifold.Contact_Count = 1

	return true
end

function Module.update_Circle_AABB(self: RigidBody)
	local position_X, position_Y = self.Position_X, self.Position_Y
	local radius = self.Vertices[1]
	self._AABB_Maximum_X, self._AABB_Maximum_Y = position_X + radius, position_Y + radius
end


-- -- -- Polygon

function Module.new_Polygon(vertices: {number}): RigidBody
	return setmetatable({
		Angle = 0,
		AngularVelocity = 0,
		Friction_Static = .5, Friction_Kinetic = .3,
		Inertia = 1, Inertia_Inverse = 1,
		Mass = 1, Mass_Inverse = 1,
		Position_X = 0, Position_Y = 0,
		Vertices = vertices,
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
local function project(verts: {number}, ax: number, ay: number): (number, number)
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
local function best_axis(a_verts: {number}, b_verts: {number}): (number, number, number, number)
	local length = #a_verts
	local best_sep = -math_Huge
	local best_nx = 0
	local best_ny = 0
	local best_edge = 1
	for index = 1, length - 2, 2 do
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

		local a_min, a_max = project(a_verts, contact_Normal_X, contact_Normal_Y)
		local b_min, b_max = project(b_verts, contact_Normal_X, contact_Normal_Y)
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
	
	
	local index = length - 1
	
	-- Compute edge normal
	local ex = a_verts[index] - a_verts[1]
	-- Safe pair indexing
	local cur_x = a_verts[index]
	local cur_y = a_verts[length]
	local next_x = a_verts[1]
	local next_y = a_verts[2]
	local edge_x = next_x - cur_x
	local edge_y = next_y - cur_y
	-- Outward normal (right-hand)
	local len = math_sqrt(edge_x * edge_x + edge_y * edge_y)
	if len < 1e-10 then
		return best_sep, best_nx, best_ny, best_edge
	end
	local contact_Normal_X = edge_y / len
	local contact_Normal_Y = -edge_x / len

	local a_min, a_max = project(a_verts, contact_Normal_X, contact_Normal_Y)
	local b_min, b_max = project(b_verts, contact_Normal_X, contact_Normal_Y)
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
	
	return best_sep, best_nx, best_ny, best_edge
end

local function findAxis(va: {number}, vb: {number}): (number, number, number, number)
	local bestSep = -math_Huge
	local bestNX, bestNY = 0, 0
	local bestEdge = 1

	for i = 1, #va, 2 do
		local j = (i + 2 <= #va) and i + 2 or 1

		local ex = va[j] - va[i]
		local ey = va[j+1] - va[i+1]

		local len = math_sqrt(ex*ex + ey*ey)
		if len < 1e-8 then continue end

		local nx = ey / len
		local ny = -ex / len

		local amin, amax = project(va, nx, ny)
		local bmin, bmax = project(vb, nx, ny)

		local sep = math_max(amin - bmax, bmin - amax)

		if sep > 0 then
			return sep, nx, ny, i -- separating axis
		end

		if sep > bestSep then
			bestSep = sep
			bestNX = nx
			bestNY = ny
			bestEdge = i
		end
	end

	return bestSep, bestNX, bestNY, bestEdge
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

-- Returns world-space vertices of a polygon
local function world_verts(self: RigidBody): {number}
	local angle = self.Angle
	local position_X, position_Y = self.Position_X, self.Position_Y
	local vertices = self.Vertices
	local length = #vertices

	local cosine, sine = math_cos(angle), math_sin(angle)

	local out: {number} = table_create(length, 0)
	for index = 1, length, 2 do
		local local_X, local_Y = vertices[index], vertices[index + 1]
		out[index], out[index + 1] = cosine * local_X - sine * local_Y + position_X, sine * local_X + cosine * local_Y + position_Y
	end

	return out
end


-- -- Methods

--[[
	Check collision between two polygons.
	returns: collided
]]
function Module.collide_Polygon(self: RigidBody, polygon_B: RigidBody, manifold: Manifold): boolean
	-- === Build world verts ===
	local va = world_verts(self)
	local vb = world_verts(polygon_B)

	-- === SAT ===
	local sepA, nxA, nyA, edgeA = findAxis(va, vb)
	if sepA > 0 then return false end

	local sepB, nxB, nyB, edgeB = findAxis(vb, va)
	if sepB > 0 then return false end

	-- === Choose reference face ===
	local refVerts, incVerts
	local refEdge
	local nx, ny

	if sepA > sepB then
		refVerts = va
		incVerts = vb
		refEdge = edgeA
		nx, ny = nxA, nyA
	else
		refVerts = vb
		incVerts = va
		refEdge = edgeB
		nx, ny = -nxB, -nyB
	end

	-- Ensure normal points A → B
	local dx = polygon_B.Position_X - self.Position_X
	local dy = polygon_B.Position_Y - self.Position_Y
	if dx*nx + dy*ny < 0 then
		nx, ny = -nx, -ny
	end

	-- === Reference edge ===
	local i1 = refEdge
	local i2 = (i1 + 2 <= #refVerts) and i1 + 2 or 1

	local v1x, v1y = refVerts[i1], refVerts[i1+1]
	local v2x, v2y = refVerts[i2], refVerts[i2+1]

	local ex = v2x - v1x
	local ey = v2y - v1y
	local len = math_sqrt(ex*ex + ey*ey)
	ex /= len
	ey /= len

	local tx, ty = ex, ey

	-- === Incident edge ===
	local bestDot = math.huge
	local incIndex = 1

	for i = 1, #incVerts, 2 do
		local j = (i + 2 <= #incVerts) and i + 2 or 1
		local ex2 = incVerts[j] - incVerts[i]
		local ey2 = incVerts[j + 1] - incVerts[i + 1]

		local l = math_sqrt(ex2*ex2 + ey2*ey2)
		if l < 1e-8 then continue end

		ex2 /= l
		ey2 /= l

		local d = ex2*nx + ey2*ny
		if d < bestDot then
			bestDot = d
			incIndex = i
		end
	end

	local j = (incIndex + 2 <= #incVerts) and incIndex + 2 or 1
	local p1x, p1y = incVerts[incIndex], incVerts[incIndex+1]
	local p2x, p2y = incVerts[j], incVerts[j+1]

	-- === Clip ===
	local offset1 = v1x * tx + v1y * ty
	local offset2 = v2x * tx + v2y * ty

	local count: number
	p1x, p1y, p2x, p2y, count = clip_segment(p1x, p1y, p2x, p2y,  tx,  ty, offset1)
	if count < 2 then return false end

	p1x, p1y, p2x, p2y, count = clip_segment(p1x, p1y, p2x, p2y, -tx, -ty, -offset2)
	if count < 2 then return false end

	-- === Final contacts ===
	local refDot = v1x * nx + v1y * ny

	local pen1 = refDot - (p1x * nx + p1y*ny)
	local pen2 = refDot - (p2x * nx + p2y*ny)

	local countFinal = 0

	if pen1 >= 0 then
		manifold.Contact1_X = p1x
		manifold.Contact1_Y = p1y
		countFinal += 1
	end

	if pen2 >= 0 then
		if countFinal == 0 then
			manifold.Contact1_X = p2x
			manifold.Contact1_Y = p2y
		else
			manifold.Contact2_X = p2x
			manifold.Contact2_Y = p2y
		end
		countFinal += 1
	end

	if countFinal == 0 then return false end

	manifold.Contact_Count = countFinal
	manifold.Contact_Normal_X = nx
	manifold.Contact_Normal_Y = ny
	manifold.Penetration = math.max(pen1, pen2)

	return true
end

function Module.update_Polygon_AABB(self: RigidBody)
	local angle = self.Angle
	local vertices = self.Vertices
	local position_X, position_Y = self.Position_X, self.Position_Y

	local cosine, sine = math_cos(angle), math_sin(angle)
	local min_x, min_y = math_Huge, math_Huge
	local max_x, max_y = -math_Huge, -math_Huge

	for index = 1, #vertices, 2 do
		local local_X, local_Y = vertices[index], vertices[index + 1]

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
