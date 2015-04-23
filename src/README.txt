HOMEWORK 1: SIMPLIFICATION & SUBDIVISION


NAME:  Austin Gulati


TOTAL TIME SPENT:  12
Please estimate the number of hours you spent on this assignment.


COLLABORATORS AND OTHER RESOURCES: List the names of everyone you
talked to about this assignment and all of the resources (books,
online reference material, etc.) you consulted in completing this
assignment.

n/a

Remember: Your implementation for this assignment must be done on your
own, as described in "Academic Integrity for Homework" handout.



GOURAUD SHADING EFFICIENCY/PERFORMANCE:

t=# of triangles
a=average # of triangles adjacent to each vertex

O(ta)

The shading algorithm requires work to be done for each triangle.
The work that is done for each triangle involves computing the normal
for each adjacent triangle. This happens three times. The runtime of
this inner step is the average number of triangles adjancent to each vertex.

SIMPLIFICATION/EDGE COLLAPSE NOTES:

I collapse the shortest edge. I find the shortest edge using the slow
linear search. It is very slow for large meshes

To check if a collapse is legal, I check the number of points in common
between the modified triangles connected to the start vertex and end vertex.
If it is >2, there is a shared vertex, and we skip that edge.
Before implementing this, I was getting some blue stuff, now I am not

When simplifying a mesh with boundaries, at some point the boundary breaks
and this messes up the whole process

SUBDIVISION NOTES:

I believe that I handle valence 6 and extraordinary vertices correctly
even though my code does nothing special to handle that. I use a B=(3/(8n))
calculation which works for both cases

I wasn't able to fully get the crease rules working

KNOWN BUGS IN YOUR CODE:
Please be concise!

Unable to simplify meshes with boundaries
Crease rules don't work


NEW FEATURES OR EXTENSIONS FOR EXTRA CREDIT:
Include instructions for use and test cases and sample output as appropriate.
