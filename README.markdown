# Enoch

## Introduction

Enoch is an object-oriented framework written in C++.  It provides an application program interface (API) that can be used within other programs to represent and analyze finite element models.  There is no user interface built into Enoch—it is up to the application programmer to gather data from users and pass it on to the Enoch code.

A structural model consists of an array of Node objects and an array of Element objects.  The Node objects store position, displacement, velocity, acceleration and mass values for an arbitrary number of degrees of freedom.  The Element class is associated with a set of Node classes and defines an interface for getting an element’s stiffness matrix.

A structural model is passed to an Enoch object that performs a Newmark-Beta Time History Analysis.  Loads are specified for each degree of freedom at each time step, allowing arbitrary, complex load conditions.  The damping matrix is generated using Rayleigh’s model of damping.

If all the elements are linear, the stiffness and damping matrices are generated at the beginning of the analysis and remain constant throughout.  Structural models incorporating non-linear elements are included in the analysis by updating the global stiffness matrix at each time step.  The linear elements are assembled into a global linear stiffness matrix, which remains constant throughout the time-history analysis.  At each time step, the tangent stiffness matrices of each nonlinear element are acquired through the interface defined in the Element class.

Non-linear elements are added by creating a new class derived from the Element class.  When an Enoch object requests the element’s stiffness matrix, any computation to generate the tangent stiffness matrix may be performed.

The abstraction of elements allows many different elements, both linear and non-linear, to be included in the same analysis.  Since the time-history analysis interacts with each element through the same interface, new elements can be developed without changing the analysis code.  Using the same interface, it is also possible to develop a new analysis module that uses the existing elements.

## Future Work

There still remains much work to be done on this program.  The data structure used to store the loading information needs to be changed.  Currently the load is stored as matrix with (number of degrees of freedom) x (number of time steps) elements.  This is very inefficient, especially for the case of a constant load.  A better solution would be to develop a Load class and associate an instance of that class with each degree of freedom.  The Load class would be inherited for different types of loads, for instance ConstantLoad, EarthquakeLoad, SinusoidalLoad, etc...  This would simplify the task of using many different loads.

The data output behaviour currently leaves much to be desired.  The ultimate solution would be for Enoch to make the data available to another class which would deal with the task of logging the data.

Finally, more validation is needed.  Once the above changes have been made, many structures should be tested in many load cases to ensure that the program gives accurate results and to gain an understanding of how small the time step needs to be to guarantee convergence.
