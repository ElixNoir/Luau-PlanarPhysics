--- Author: COMT

--[[

	PlanarPhysics - A module for handling planar physics.
	Mod: Manifold
	Version: 0.0.1

]]


-- -- State

local Module = true

export type Manifold = {
	Bias: number,
	Contact_Count: number,
	Contact_Normal_X: number,
	Contact_Normal_Y: number,
	Contact1_X: number, Contact1_Y: number,
	Contact2_X: number, Contact2_Y: number,
	Impulse_n1: number, Impulse_t1: number,
	Impulse_n2: number, Impulse_t2: number,
	Mass_n: number, Mass_t: number,
	Penetration: number,
	RigidBodyIdentity_A: number,
	RigidBodyIdentity_B: number,
	_Fresh: boolean,
	_Key: number
}


-- -- End

return Module :: typeof(Module) & {
	Type: Manifold
}
