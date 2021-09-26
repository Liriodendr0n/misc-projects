This script takes a non-traditional approach to aircraft sizing, leveraging the 'enormous' computing power we have access to in modern personal computers.
It takes a multidimensional parameter space, and evaluates the performance of every single element of that parameter space.
Each configuration is then checked against the relevant constraints (TOFL, battery energy, ...), and the configurations with the highest scores are deemed the winners.
My laptop can process around 30,000 configuration candidates per second, your performance may vary.

Do not make critical decisions based on this code's results yet, it still needs vastly improved parameter estimation methods.

Methods Wanted:

DRAG CALCULATION
oswald efficiency as a function of wing geometry, airfoil, and streamlining
CLmax as a function of airfoil and wing geometry

WEIGHT ESTIMATION
anything, anything better than a multiple of the payload weight (right now it's just some factor times payload weight)


PLOT INFORMATION
Right now contour plots are made with the mission 2 and 3 scores as the X and Y axes, and a design variable as the height/color
The upper-rightmost border of the plot is the pareto front, the locus of configurations where no individual score can be improved without reducing another. 
Those are the configurations we want to decide from, there may be many, there may be one, they will depend on constraints, ranges, and modeling assumptions.
Visual analysis of the plots can determine which design variables are "independent" and "dependent".
"Independent" design variables are the "dials" we can turn that affect our score in comprehensible ways, whereas "dependent" variables are "noisy dials" that do not have easily comprehensible effects on the score, and are determined by constraints and independent variables.