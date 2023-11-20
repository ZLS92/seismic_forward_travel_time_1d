# seismic_forward_travel_time_1d
Python code to compute the travel time and the raypath of seismic wases direct arrivals, crossing a 1d media.
The code has been implemented following the solution described in the user manual of the softwaer Cat3d by Gualtiero Böhm (see refernces).

## Installation
The code requires only the installation of the `numpy` and `matplotlib` libraries. 
Any version of `numpy` and `matplotlib` compatible with Python 3 should be sufficient.
To run the code, simply execute it from your terminal, for example:
```
> python seismic_forward_travel_time_1D.py
```

## Minimum time path
This algorithm is designed to calculate the minimum time path between a source point and a receiver point within a 1D velocity model. 
It utilizes an iterative procedure that operates on individual triplets of points, representing intersections of the ray with the velocity layers. 
The process begins with the initial assumption that, in the case of a directly transmitted wave, a straight-line segment connects the source to the receiver.

The algorithm unfolds in three iterative steps:

1. Iteration to determine the minimum travel time for each triplet of raypath points (source - interface intersection - receiver).
2. Iteration over each triplet crossing a specified interface at depth.
3. Iteration over each specified source-receiver pair.

### First step, i.e., "def min_travel_time_SinglePlane(...)" 

```
   R
   |                                         x = |Rp-I|                                 
   ||                                        a = |R-Rp|                                       
   | |                                       b = |S-Sp|                                       
   |  |                                      c = |Sp-Rp|
   |   |                                    
   |    |                                   
   |     |                  V1
 a |      |                                                                         
   |       |                                                                        
   |        |                                                                       
   |         |                                                                      
   |          |                   
   |           |                  
   |            |                  
 Rp|             | I                c                          plane p
----------------------------------------------------------------------
         x         |                                   |Sp
                      |                                |                        
                         |                             |                     
                            |                          |                   
                               |                       |                   
                                  |                    |                    
                                     |                 | b
                                        |              |       
                                           |           |             
                             V2               |        |                
                                                 |     |                   
                                                    |  |
                                                       |
                                                       S
```

Given the triplet R-I-S, with R and S located in two half-spaces separated by the plane p, the problem is to find the path from S to R in the shortest possible time, 
taking into account the velocities V1 and V2 in the two half-spaces.
After simple calculations (using the Pythagorean theorem), the time for the R-I-S path is given by the function:

$$ F(x) = \frac{\sqrt{a^2 + x^2}}{V_1} + \frac{\sqrt{b^2 + (c - x)^2}}{V_2} $$

To find the minimum of this function, we examine where its first derivative equals zero $( F'(x)=0 )$.

$$ F'(x) = \frac{x}{V_1 \sqrt{a^2 + x^2}} - \frac{c - x}{V_2 \sqrt{b^2 + (c - x)^2}} $$

At this point, the bisection method is used to find the value of $x$ for which the function becomes zero. The values of the derivative $F1$​, $F2$​, and $F3​$ 
are calculated for the three points (x=0, x=c/2, and x=c). 
It is checked where the zero lies between  $F1$​, $F2$​, and $F3​$; the excluded value is discarded, and a new derivative value is calculated by halving x.
For example, if $F1>0$​, $F2>0$​, and $F3<0​$, $F1$​ is discarded, and a new $F4$​ is calculated for $x=2/3​c$. 
The same reasoning is applied to the new values $F2$​, $F3$​, and $F4$​, and so on until the value $Fn$​ is less than a pre-set threshold (very small, of course).

### Second step, i.e., "def min_travel_time_MultiPlane(...)"

The second step iterate over each triplet of points of the raypath. For example, the second triplet considers the previusly found intersection as source point (the new S is the I from previus step). 
This process continues until reaching the receiver and then retracing the path back to the source with the same procedure, 
repeating these steps until the difference between the radius of the previous step and the last calculated radius is less than a certain threshold. 
This iterative procedure is applied to all intersections of the ray segment.

### Third step, i.e., "def def min_travel_times_paths_1D(...)"
The last step iterate the second step over each couple of the input source-reciever points, returning a list of all the corresponding travel-times and ray-paths.
Optionally, it is possible to specify which couple of source and reciver input to compute by providing the name of the sources and of the recivers in the fourth column of the source_poits, reciever_points arrays, 
and also providing the aegument SR_cobinations with a 2 columns array with the name of the selected sources and recievers on each line.

### References
* *Böhm, G. and OGS research group (2014), Cat3D - Computer Aided Tomography for 3-D models. User Manual, OGS; https://hdl.handle.net/20.500.14083/6558*
* *Böhm, G., G.Rossi, e A. Vesnaver, 1999. Minimum time ray-tracing for 3-D irregular grids, J. of Seism. Expl., 8: 117-131.*




