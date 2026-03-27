--- Author: COMT

--[[

	PlanarPhysics - A module for handling planar physics.
	Mod: World
	Version: 0.0.1

]]


--!native
--!optimize 2
--!strict


-- -- State

local Manifold = require("./Manifold")
type Manifold = Manifold.Manifold

local RigidBody = require("./RigidBody")
type RigidBody = RigidBody.RigidBody

local Module = {}
Module.__index = Module

type Island = {
	Manifolds: {number},
	RigidBodyIdentities: {number},
	Sleeping: boolean
}

export type World = typeof(setmetatable({} :: {
	Acceleration_X: number,
	Acceleration_Y: number,
	Baumgarte: number,
	_BroadPhase_Pairs: {number},
	_BroadPhase_Sorted: {number},
	Iterations: number,
	RigidBodies: {RigidBody},
	Sleep_Threshold_V: number,
	Sleep_Threshold_W: number,
	Sleep_Time_Required: number,
	Slop: number,
	_Islands: {Island},
	_Manifolds: {Manifold},
	_Manifold_Map: {[number]: Manifold}
}, Module))


-- -- Helpers

local math_abs = math.abs
local math_clamp = math.clamp
local math_max = math.max
local math_min = math.min
local math_sqrt = math.sqrt

local table_clear = table.clear
local table_create = table.create
local table_insert = table.insert
local table_remove = table.remove
local table_sort = table.sort

local Pi2 = 2 * math.pi

local SLEEP_V_SQ = 0.0001 -- squared linear speed threshold
local SLEEP_W_SQ = 0.0001 -- squared angular speed threshold
local SLEEP_TIME = 0.5 -- seconds at rest before sleeping
local BAUMGARTE = 0.2
local SLOP = 0.005 -- penetration slop (meters)
local MAX_ITERS = 10

local RigidBody_collide_Circle_Circle = RigidBody.collide_Circle_Circle
local RigidBody_collide_Circle_OBB = RigidBody.collide_Circle_OBB
local RigidBody_collide_Polygon = RigidBody.collide_Polygon

local function broadphase_pairs(self: World)
	local broadPhase_Pairs = self._BroadPhase_Pairs
	local broadPhase_Sorted = self._BroadPhase_Sorted
	local rigidBodies = self.RigidBodies
	local rigidBodies_Length = #rigidBodies

	table_clear(broadPhase_Sorted)
	for index = 1, rigidBodies_Length do
		broadPhase_Sorted[index] = index
	end
	
	table_sort(broadPhase_Sorted, function(index_A: number, index_B: number): boolean
		return rigidBodies[index_A]._AABB_Minimum_X < rigidBodies[index_B]._AABB_Minimum_X
	end)

	local pair_Count = 1
	for i = 1, rigidBodies_Length do
		local rigidBody_A_Identity = broadPhase_Sorted[i]
		local rigidBody_A = rigidBodies[rigidBody_A_Identity]
		if rigidBody_A._Sleeping then
			continue
		end
		
		for j = i + 1, rigidBodies_Length do
			local rigidBody_B_Identity = broadPhase_Sorted[j]
			local rigidBody_B = rigidBodies[rigidBody_B_Identity]
			if rigidBody_B._AABB_Minimum_X > rigidBody_A._AABB_Maximum_X then
				break
			end
			
			if (rigidBody_B._Sleeping and rigidBody_A.Mass_Inverse == 0)
				or (rigidBody_A._AABB_Maximum_Y < rigidBody_B._AABB_Minimum_Y)
				or (rigidBody_A._AABB_Minimum_Y > rigidBody_B._AABB_Maximum_Y)
				or (rigidBody_A.Mass_Inverse == 0 and rigidBody_B.Mass_Inverse == 0) then
				continue
			end
			
			broadPhase_Pairs[pair_Count], broadPhase_Pairs[pair_Count + 1] = rigidBody_A_Identity, rigidBody_B_Identity
			pair_Count += 2
		end
	end
	
	for index = #broadPhase_Pairs, pair_Count + 1, -1 do
		broadPhase_Pairs[index] = nil
	end
end

local function compute_manifold(rigidBody_A: RigidBody, rigidBody_B: RigidBody, manifold: Manifold): boolean
	local ca, cb = #rigidBody_A.Shape == 1, #rigidBody_B.Shape == 1

	if ca and cb then
		return RigidBody_collide_Circle_Circle(rigidBody_A, rigidBody_B, manifold)
	elseif ca then
		return RigidBody_collide_Circle_OBB(rigidBody_A, rigidBody_B, manifold, true)
	elseif cb then
		return RigidBody_collide_Circle_OBB(rigidBody_B, rigidBody_A, manifold, false)
	end
	
	return RigidBody_collide_Polygon(rigidBody_A, rigidBody_B, manifold)
end

local function contact_manager_update(self: World)
	local broadPhase_Pairs = self._BroadPhase_Pairs
	local rigidBodies = self.RigidBodies

	local oldMap = self._Manifold_Map

	local newManifolds: {Manifold} = {}
	local newMap: {[number]: Manifold} = {}

	local manifoldCount = 0

	for index = 1, #broadPhase_Pairs, 2 do
		local rigidBodyIdentity_A = broadPhase_Pairs[index]
		local rigidBodyIdentity_B = broadPhase_Pairs[index + 1]

		local key = rigidBodyIdentity_A <= rigidBodyIdentity_B and rigidBodyIdentity_A * 2097152 + rigidBodyIdentity_B or rigidBodyIdentity_B * 2097152 + rigidBodyIdentity_A

		local manifold = oldMap[key]

		if not manifold then
			manifold = {
				RigidBodyIdentity_A = rigidBodyIdentity_A,
				RigidBodyIdentity_B = rigidBodyIdentity_B,

				Contact_Normal_X = 0,
				Contact_Normal_Y = 0,

				Penetration = 0,

				Contact1_X = 0,
				Contact1_Y = 0,
				Contact2_X = 0,
				Contact2_Y = 0,

				Contact_Count = 0,

				Impulse_n1 = 0,
				Impulse_t1 = 0,
				Impulse_n2 = 0,
				Impulse_t2 = 0,

				Mass_n = 0,
				Mass_t = 0,
				Bias = 0,

				_Fresh = false,
				_Key = key
			}
		end

		if compute_manifold(rigidBodies[rigidBodyIdentity_A], rigidBodies[rigidBodyIdentity_B], manifold) then
			local rigidBody_A = rigidBodies[rigidBodyIdentity_A]
			local rigidBody_B = rigidBodies[rigidBodyIdentity_B]

			if rigidBody_A._Sleeping then
				rigidBody_A._Sleeping = false
				rigidBody_A._Sleep_Timer = 0
			end

			if rigidBody_B._Sleeping then
				rigidBody_B._Sleeping = false
				rigidBody_B._Sleep_Timer = 0
			end

			manifoldCount += 1
			newManifolds[manifoldCount] = manifold
			newMap[key] = manifold
		end
	end

	for i = manifoldCount + 1, #newManifolds do
		newManifolds[i] = nil
	end

	self._Manifolds = newManifolds
	self._Manifold_Map = newMap
end

local function integrate(self: World, deltaTime: number)
	local acceleration_X, acceleration_Y = self.Acceleration_X, self.Acceleration_Y
	for _, rigidBody in self.RigidBodies do
		if rigidBody._Sleeping or rigidBody.Mass_Inverse == 0 then
			continue
		end
		
		-- Semi-implicit Euler
		if rigidBody.Angle then
			rigidBody.Angle = (rigidBody.Angle + rigidBody.AngularVelocity * deltaTime) % Pi2
		end
		
		local velocity_X, velocity_Y = rigidBody.Velocity_X + acceleration_X * deltaTime, rigidBody.Velocity_Y + acceleration_Y * deltaTime
		rigidBody.Position_X += velocity_X * deltaTime
		rigidBody.Position_Y += velocity_Y * deltaTime
		rigidBody.Velocity_X = velocity_X
		rigidBody.Velocity_Y = velocity_Y
		rigidBody:update_AABB()
	end
end

local function rebuild_Islands(self: World)
	local rigidBodies = self.RigidBodies
	local manifolds = self._Manifolds
	local rigidBodies_Length = #rigidBodies

	local parent: {number} = table_create(rigidBodies_Length, 0)
	for index = 1, rigidBodies_Length do
		parent[index] = index
	end

	for _, manifold in manifolds do
		local rx = manifold.RigidBodyIdentity_A
		while parent[rx] ~= rx do
			parent[rx] = parent[parent[rx]]
			rx = parent[rx]
		end
		
		local ry = manifold.RigidBodyIdentity_B
		while parent[ry] ~= ry do
			parent[ry] = parent[parent[ry]]
			ry = parent[ry]
		end
		
		if rx ~= ry then
			parent[rx] = ry
		end
	end

	local island_map: {[number]: Island} = {}
	for index = 1, rigidBodies_Length do
		local root = index
		while parent[root] ~= root do
			parent[root] = parent[parent[root]]
			root = parent[root]
		end
		
		local rigidBodyIdentities
		if not island_map[root] then
			rigidBodyIdentities = {}
			island_map[root] = {Manifolds = {}, RigidBodyIdentities = rigidBodyIdentities, Sleeping = false}
		else
			rigidBodyIdentities = island_map[root].RigidBodyIdentities
		end
		
		table_insert(rigidBodyIdentities, index)
	end
	
	for manifoldIdentity, manifold in manifolds do
		local root = manifold.RigidBodyIdentity_A
		while parent[root] ~= root do
			parent[root] = parent[parent[root]]
			root = parent[root]
		end
		
		local island = island_map[root]
		if island then
			table_insert(island.Manifolds, manifoldIdentity)
		end
	end

	local islands = {}
	self._Islands = islands
	for _, island in island_map do
		table_insert(islands, island)
	end
end

local function solver_iterate(self: World)
	local rigidBodies = self.RigidBodies
	local manifolds = self._Manifolds

	for _, manifold in manifolds do
		local rigidBody_A, rigidBody_B = rigidBodies[manifold.RigidBodyIdentity_A], rigidBodies[manifold.RigidBodyIdentity_B]
		local contact_Normal_X, contact_Normal_Y = manifold.Contact_Normal_X, manifold.Contact_Normal_Y
		local manifold_Bias = manifold.Bias
		local manifold_Mass_n = manifold.Mass_n
		
		local rigidBody_A_AngularVelocity = rigidBody_A.AngularVelocity
		local rigidBody_A_Friction_Kinetic = rigidBody_A.Friction_Kinetic
		local rigidBody_A_Inertia_Inverse = rigidBody_A.Inertia_Inverse
		local rigidBody_A_Mass_Inverse = rigidBody_A.Mass_Inverse
		local rigidBody_A_Position_X = rigidBody_A.Position_X
		local rigidBody_A_Position_Y = rigidBody_A.Position_Y
		local rigidBody_A_Velocity_X = rigidBody_A.Velocity_X
		local rigidBody_A_Velocity_Y = rigidBody_A.Velocity_Y
		
		local rigidBody_B_AngularVelocity = rigidBody_B.AngularVelocity
		local rigidBody_B_Friction_Kinetic = rigidBody_B.Friction_Kinetic
		local rigidBody_B_Inertia_Inverse = rigidBody_B.Inertia_Inverse
		local rigidBody_B_Mass_Inverse = rigidBody_B.Mass_Inverse
		local rigidBody_B_Position_X = rigidBody_B.Position_X
		local rigidBody_B_Position_Y = rigidBody_B.Position_Y
		local rigidBody_B_Velocity_X = rigidBody_B.Velocity_X
		local rigidBody_B_Velocity_Y = rigidBody_B.Velocity_Y

		local function resolve_point(contact_X: number, contact_Y: number, accumulatedNormalImpulse: number, accumulatedTangentImpulse: number): (number, number)
			-- Relative positions from centers to contact

			local raX = contact_X - rigidBody_A_Position_X
			local raY = contact_Y - rigidBody_A_Position_Y

			local rbX = contact_X - rigidBody_B_Position_X
			local rbY = contact_Y - rigidBody_B_Position_Y


			-- Relative velocity at contact point

			local relativeVelocityX =
				rigidBody_B_Velocity_X - rigidBody_A_Velocity_X
				+ rigidBody_A_AngularVelocity * raY - rigidBody_B_AngularVelocity * rbY

			local relativeVelocityY =
				rigidBody_B_Velocity_Y - rigidBody_A_Velocity_Y
				+ rigidBody_B_AngularVelocity * rbX - rigidBody_A_AngularVelocity * raX

			-- Normal direction (collision axis)

			local velocityNormal = relativeVelocityX * contact_Normal_X + relativeVelocityY * contact_Normal_Y


			-- Compute normal impulse

			local deltaNormalImpulse = manifold_Mass_n * (-velocityNormal + manifold_Bias)

			local newNormalImpulse = math_max(accumulatedNormalImpulse + deltaNormalImpulse, 0)

			deltaNormalImpulse = newNormalImpulse - accumulatedNormalImpulse


			-- Tangent direction (friction axis)

			local tangentX = -contact_Normal_Y
			local tangentY = contact_Normal_X
			local velocityTangent = relativeVelocityX * tangentX + relativeVelocityY * tangentY


			-- Compute friction impulse

			local deltaTangentImpulse = manifold_Mass_n * (-velocityTangent)

			local frictionCoefficient =
				(rigidBody_A_Friction_Kinetic + rigidBody_B_Friction_Kinetic) * 0.5

			local maxFrictionImpulse = frictionCoefficient * newNormalImpulse

			local newTangentImpulse = math_clamp(
				accumulatedTangentImpulse + deltaTangentImpulse,
				-maxFrictionImpulse,
				maxFrictionImpulse
			)

			deltaTangentImpulse = newTangentImpulse - accumulatedTangentImpulse


			-- Combine impulses into world-space vector

			local impulseX =
				deltaNormalImpulse * contact_Normal_X +
				deltaTangentImpulse * tangentX

			local impulseY =
				deltaNormalImpulse * contact_Normal_Y +
				deltaTangentImpulse * tangentY


			-- Apply impulse to both rigid bodies

			rigidBody_A.Velocity_X -= impulseX * rigidBody_A_Mass_Inverse
			rigidBody_A.Velocity_Y -= impulseY * rigidBody_A_Mass_Inverse
			rigidBody_A.AngularVelocity -= (raX * impulseY - raY * impulseX) * rigidBody_A_Inertia_Inverse

			rigidBody_B.Velocity_X += impulseX * rigidBody_B_Mass_Inverse
			rigidBody_B.Velocity_Y += impulseY * rigidBody_B_Mass_Inverse
			rigidBody_B.AngularVelocity += (rbX * impulseY - rbY * impulseX) * rigidBody_B_Inertia_Inverse

			return newNormalImpulse, newTangentImpulse -- Return updated accumulated impulses (for warm starting)
		end

		manifold.Impulse_n1, manifold.Impulse_t1 = resolve_point(manifold.Contact1_X, manifold.Contact1_Y, manifold.Impulse_n1, manifold.Impulse_t1)
		if manifold.Contact_Count == 2 then
			manifold.Impulse_n2, manifold.Impulse_t2 = resolve_point(manifold.Contact2_X, manifold.Contact2_Y, manifold.Impulse_n2, manifold.Impulse_t2)
		end
	end
end

local function solver_pre_step(self: World, inv_dt: number)
	local rigidBodies = self.RigidBodies
	local manifolds = self._Manifolds
	local baumgarte = self.Baumgarte
	local slop = self.Slop

	for _, manifold in manifolds do
		local rigidBody_A, rigidBody_B = rigidBodies[manifold.RigidBodyIdentity_A], rigidBodies[manifold.RigidBodyIdentity_B]
		local contact_Normal_X, contact_Normal_Y = manifold.Contact_Normal_X, manifold.Contact_Normal_Y

		local function point_mass(px: number, py: number): number
			local rax, ray = px - rigidBody_A.Position_X, py - rigidBody_A.Position_Y
			local rbx, rby = px - rigidBody_B.Position_X, py - rigidBody_B.Position_Y
			local ra_cross_n, rb_cross_n = rax * contact_Normal_Y - ray * contact_Normal_X, rbx * contact_Normal_Y - rby * contact_Normal_X
			local denominator = rigidBody_A.Mass_Inverse + rigidBody_B.Mass_Inverse + rigidBody_A.Inertia_Inverse * ra_cross_n * ra_cross_n + rigidBody_B.Inertia_Inverse * rb_cross_n * rb_cross_n
			return denominator > 1e-10 and 1 / denominator or 0
		end

		manifold.Mass_n = point_mass(manifold.Contact1_X, manifold.Contact1_Y)
		manifold.Bias = -baumgarte * inv_dt * math_min(0, -manifold.Penetration + slop)

		local function apply_impulse(px: number, py: number, normal_Impulse: number, jt: number)
			local normal_Impulse_X = normal_Impulse * contact_Normal_X - jt * contact_Normal_Y
			local normal_Impulse_Y = normal_Impulse * contact_Normal_Y + jt * contact_Normal_X
			local rax, ray = px - rigidBody_A.Position_X, py - rigidBody_A.Position_Y
			local rbx, rby = px - rigidBody_B.Position_X, py - rigidBody_B.Position_Y
			rigidBody_A.Velocity_X -= normal_Impulse_X * rigidBody_A.Mass_Inverse
			rigidBody_A.Velocity_Y -= normal_Impulse_Y * rigidBody_A.Mass_Inverse
			rigidBody_A.AngularVelocity -= (rax * normal_Impulse_Y - ray * normal_Impulse_X) * rigidBody_A.Inertia_Inverse
			rigidBody_B.Velocity_X += normal_Impulse_X * rigidBody_B.Mass_Inverse
			rigidBody_B.Velocity_Y += normal_Impulse_Y * rigidBody_B.Mass_Inverse
			rigidBody_B.AngularVelocity += (rbx * normal_Impulse_Y - rby * normal_Impulse_X) * rigidBody_B.Inertia_Inverse
		end

		apply_impulse(manifold.Contact1_X, manifold.Contact1_Y, manifold.Impulse_n1, manifold.Impulse_t1)
		if manifold.Contact_Count == 2 then
			apply_impulse(manifold.Contact2_X, manifold.Contact2_Y, manifold.Impulse_n2, manifold.Impulse_t2)
		end
	end
end

local function update_sleep(self: World, deltaTime: number)
	local rigidBodies = self.RigidBodies
	for _, island in self._Islands do
		local can_sleep = true
		for _, rigidBodyIdentity in island.RigidBodyIdentities do
			local rigidBody = rigidBodies[rigidBodyIdentity]
			if rigidBody.Mass_Inverse == 0 then
				continue
			end
			local b_AngularVelocity = rigidBody.AngularVelocity
			local b_Velocity_X, b_Velocity_Y = rigidBody.Velocity_X, rigidBody.Velocity_Y
			if b_Velocity_X * b_Velocity_X + b_Velocity_Y * b_Velocity_Y > SLEEP_V_SQ or b_AngularVelocity * b_AngularVelocity > SLEEP_W_SQ then
				rigidBody._Sleep_Timer = 0
				can_sleep = false
			else
				rigidBody._Sleep_Timer += deltaTime
				if rigidBody._Sleep_Timer < SLEEP_TIME then can_sleep = false end
			end
		end
		if can_sleep then
			island.Sleeping = true
			for _, bi in island.RigidBodyIdentities do
				local b = rigidBodies[bi]
				b._Sleeping = true
				b.Velocity_X, b.Velocity_Y = 0, 0
				b.AngularVelocity = 0
			end
		end
	end
end


-- -- Functions

function Module.new(): World
	return setmetatable({
		Acceleration_X = 0,
		Acceleration_Y = 9.81,
		Baumgarte = BAUMGARTE,
		Iterations = MAX_ITERS,
		RigidBodies = {},
		Sleep_Threshold_V = math_sqrt(SLEEP_V_SQ),
		Sleep_Threshold_W = math_sqrt(SLEEP_W_SQ),
		Sleep_Time_Required = SLEEP_TIME,
		Slop = SLOP,
		_BroadPhase_Pairs = {},
		_BroadPhase_Sorted = {},
		_Islands = {},
		_Manifolds = {},
		_Manifold_Map = {}
	}, Module)
end


-- -- Methods

function Module.addRigidBody(self: World, rigidBody: RigidBody)
	local rigidBodies = self.RigidBodies
	local index = #rigidBodies + 1
	rigidBody._Index = index
	rigidBodies[index] = rigidBody
end

function Module.removeRigidBody(self: World, rigidBody: RigidBody)
	local index = rigidBody._Index
	if index == 0 then
		return
	end
	
	local rigidBodies = self.RigidBodies
	local last = #rigidBodies
	if index ~= last then
		rigidBodies[index] = rigidBodies[last]
		rigidBodies[index]._Index = index
	end
	
	rigidBodies[last] = nil
	rigidBody._Index = 0
	
	local manifolds = self._Manifolds
	local manifold_Map = self._Manifold_Map
	for i = #manifolds, 1, -1 do
		local manifold = manifolds[i]
		local rigidBody_A, rigidBody_B = manifold.RigidBodyIdentity_A, manifold.RigidBodyIdentity_B
		if rigidBody_A == last or rigidBody_B == last or rigidBody_A == index or rigidBody_B == index then
			manifold_Map[manifold._Key] = nil
			table_remove(manifolds, i)
		end
	end
end

function Module.step(self: World, deltaTime: number)
	if deltaTime <= 0 then
		return
	end

	local inv_dt = 1 / deltaTime

	-- 1. Integrate forces → update velocities & positions
	integrate(self, deltaTime)

	-- 2. Broadphase
	broadphase_pairs(self)

	-- 3. Narrowphase + contact manager
	contact_manager_update(self)

	-- 4. Island / sleep bookkeeping
	rebuild_Islands(self)

	-- 5. Pre-step (warm start + constraint mass)
	solver_pre_step(self, inv_dt)

	-- 6. Iterative impulse solver
	for _ = 1, self.Iterations do
		solver_iterate(self)
	end

	-- 7. Sleep update
	update_sleep(self, deltaTime)
end


-- -- End

return Module
