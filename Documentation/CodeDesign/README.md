## This is a todo list

- ~~Create Cell Based derived Quantities~~
- ~~Link Derived Quantities to Recruitment Processes~~
- Figure out how to set scalars on agents when linking modeled number of agents to weight from the collective individuals
- Currently derived quantities wrap a mortality block following Casal2 formulations, I am thinking of pulling it out and creating it as a normal process? 
 This means that recruitment by definition has to follow a mortality block (because derived quantities are associated with mortality blocks), which seems like a shit limitation.
- Spatial mortality
- Spatial movement
- Currently WorldCell sets agents, but that means I have to give it access to processes or put stuff on the model so it can access. I chose this route because
it means that only WorldCell has access to the Agent class (Some sort of encapsulation), I didn't want to give all the processes and initialisation phases
the ability to create agents. I am starting to think this is too restrictive though, in terms of expanding agent definition