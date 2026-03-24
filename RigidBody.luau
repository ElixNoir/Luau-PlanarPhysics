--- Author: COMT

--[[

  PlanarPhysics - A module for handling planar physics.
  Mod: RigidBody

]]

--!native
--!optimize 2
--!strict


-- -- State

local Module = {}
Module.__index = Module

export type RigidBody = typeof(setmetatable({} :: {
  Angle: number,
  Angle_Velocity: number,
  Inertia: number,
  Inverted_Inertia: number,
  Inverted_Mass: number,
  Mass: number,
  Position_X: number,
  Position_Y: number,
  Shape: {number},
  Velocity_X: number,
  Velocity_Y: number
}, Module))


-- -- Functions

function Module.new(): RigidBody
  return setmetatable({
    Angle = 0,
    Angle_Velocity = 0,
    Inertia = 1,
    Inverted_Inertia = 1,
    Inverted_Mass = 1,
    Mass = 1,
    Position_X = 0,
    Position_Y = 0,
    Shape = {},
    Velocity_X = 0,
    Velocity_Y = 0
  }, Module)
end


-- -- Methods

function Module.integrate(self: RigidBody, deltaTime: number)
  self.Position_X += self.Velocity_X * deltaTime
  self.Position_Y += self.Velocity_Y * deltaTime
  self.Angle += self.Angle_Velocity * deltaTime
end


-- -- End

return Module
