![WetPlanets_by_RoberOost_and_AnnemarieJutte](https://user-images.githubusercontent.com/50207671/114404412-f78c2c80-9ba5-11eb-82ee-77d0e7fb0047.png)

## Game Physics - Mini Project
### 12 April 2021


This repository contains the Unity project created for the mini project of the course Game Physics. In our project we created a simulation modeling fluid flow as a velocity field on a MAC grid. The velocity field was modeled using the Navier-Stokes equations. The implementation is based on the paper Fluid Flow for the Rest of Us by Cline et al [1]. In the simulation the fluid is placed on a planet and the influence of an orbiting planet is shown. The influence of the planets on the fluid is found using Newton's law, as is the planet's orbit itself. The planetary physics are based on Lague's tutorial [2]. 


## Instructions

The project can be run in Unity, the scene PlanetFluidScene contains the main simulation. The scene FluidScene contains a basic fluid simulation, without the planetary influences and the scene PlanetScene contains a simulation of only planetary orbits. When running the PlanetFluidScene in Unity various parameters can be tweaked in the inspector by selecting the right game objects. The SimulationMaster object contains the most interesting parameters influencing the fluid simulation.


## References
[1] Cline, D., Cardon, D., & Egbert, P. K. (2013). Fluid flow for the rest of us: Tutorial of the marker and cell method in computer graphics. Brigham Young University.

[2] Lague, S., Coding Adventure: Solar System (2020). https://github.com/SebLague/Solar-System/tree/Episode_01

## Assets Used
- **Planet materials:** Vast Outer Space by Prodigious Creations. https://assetstore.unity.com/packages/3d/environments/sci-fi/vast-outer-space-38913
- **Skybox:** MilkyWay by Adam Bielecki https://assetstore.unity.com/packages/2d/textures-materials/milky-way-skybox-94001
