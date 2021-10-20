## Requirements 
You need to have minizinc installed. 
You need to have Julia installed to compiled the solveur. 

Include SeaPearl in minizinc solver

### 1 Building the project
Run in the console 
`julia compiler.jl`

It might take a while.
It will create a new folder named Compiled. 

### 2 Add SeaPearl to Minizinc IDE
Open minizinc IDE.
Then go to Minizinc IDE -> Preferences
In the solver box, click on add new.
Give the name, solver ID and version you want. 
In the executable part, add this this folder : 
'<path_to_the_dossier>/Compiled/bin/Exec'
In the library path, add this folder : 
'<path_to_the_dossier>/Library'
Click on add. 

