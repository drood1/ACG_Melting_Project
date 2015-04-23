HOMEWORK 2: CLOTH & FLUID SIMULATION

NAME: Austin Gulati


TOTAL TIME SPENT: 30
Please estimate the number of hours you spent on this assignment.


COLLABORATORS AND OTHER RESOURCES: List the names of everyone you
talked to about this assignment and all of the resources (books,
online reference material, etc.) you consulted in completing this
assignment.

Remember: Your implementation for this assignment must be done on your
own, as described in "Academic Integrity for Homework" handout.


DISCUSS STABILITY OF YOUR CLOTH SIMULATION:

Run with timesteps 10x:
./simulation -cloth ../src/small_cloth.txt -timestep 0.01
./simulation -cloth ../src/provot_original.txt -timestep 0.01
./simulation -cloth ../src/provot_correct_structural.txt -timestep 0.01
./simulation -cloth ../src/provot_correct_structural_and_shear.txt -timestep 0.01

The velocities seem to get out of control and the cloth wiggles around a bit but other than that it's fairly stable.

DESCRIBE YOUR NEW CLOTH TEST SCENE & 
THE COMMAND LINE TO RUN YOUR EXAMPLE:

It is similar to the square table but is a round table

./simulation -cloth ../src/table_round.txt -timestep 0.01

DISCUSS THE ACCURACY & STABILITY OF YOUR FLUID SIMULATION:

Fairly accurate, one water point doesn't move in the drop simulation but other than that I think things are working normally.
Velocities get out of control eventually on some of the simulations.

KNOWN BUGS IN YOUR CODE:
Please be concise!

Timesteps adjusted on cloth
One water point does not move on fluid_drop.txt


NEW FEATURES OR EXTENSIONS FOR EXTRA CREDIT:
Include instructions for use and test cases and sample output as appropriate.

N/A