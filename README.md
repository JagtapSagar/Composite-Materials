# Composite-Materials-Project

A composite laminate contains any number of composite lamina arranged at diferent angles.
A lamina is typically modeled as having one layer of fibers through the thickness.

The aim is to find the optimum thickness of each layer/lamina of the laminate in such a way that none of the layers fail. For each layer we consider top, middle, and bottom planes to calculate the Strength Ratio and the Tsai Wu value.

Counter intuitive to minimize the failure index when designing a laminate under this constraint; where results would be load- dependent. Maximization of the safety factor instead yields load-independent results.

In a uni-axial loading situation, predicting failure obviously reduces itself to comparing the internal stresses (s) to the material’s strength (F) in the loading direction. In this situation, the
failure index (FI) and the strength ratio (SR) are defined as:

FI = sigma/F

SR = F/sigma

In this simple scenario, the failure and strength ratio are clearly related as they are the inverse of one another. Extending this example to a composite material layer and comparing each component of the stress state to the material strength properties, one can draw the ply failure envelope for the maximum stress criterion.
