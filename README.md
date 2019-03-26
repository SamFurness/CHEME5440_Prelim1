# CHEME5440_Prelim1
Sam Furness
3/26/19


Question1:
  All details are in SamFurness_CHEME5440Prelim#1.pdf

Question2:

  Part A is solved with CHEME5440P1_2a.m and gives the plot. It calls the function bind_fxn.m to solve each fI.
  
  Part B is solved in CHEME5440P1_2bc.m. This calls the functions compute.m and bind_fxn.m. compute.m is basically Part A transformed into   a function. It does not print the sensitivity matrices (too much data). If you wish to look at a sensitivity array, the variable names     are "s_'name of variable for that parameter'" This gives you the sensitivity coefficient at each time for that parameter (note that       data is only in the phases analyzed).
  
  Part C is solved in CHEME5440P1_2bc.m and the first columns of the U matrix are returned for each phase.
  
  NOTE: SamFurness_CHEME5440Prelim#1.pdf has more information on all methods and answers to the questions
  
Question 3:

  Part A is explained in SamFurness_CHEME5440Prelim#1.pdf and the stoichiometric matrix is in CHEME5440P1_3a,b.jl
  
  Part B is solved with CHEME5440P1_3a,b.jl. It calls the function Flux.jl.
  
  Part C is solved with CHEME5440P1_3c.jl. It calls the function Flux.jl. It prints out shadow prices analyzed in                           SamFurness_CHEME5440Prelim#1.pdf
  
  NOTE: SamFurness_CHEME5440Prelim#1.pdf has more information on all methods and answers to the questions
  
  
  
