include("parser.jl")
using SeaPearl
    using Pkg

mutable struct Interpreter
    node::AST
    GLOBAL_VARIABLE::Dict
    GLOBAL_CONSTRAINT::Array
    GLOBAL_PARAMETER::Dict
    MINIMIZATION::Bool
    output::String
    Interpreter(node) = new(node, Dict(), [], Dict(), true, "")
end

function error(message)
    throw(ArgumentError(message))
end


function create_variable(interpreter::Interpreter, node::AST, trailer, m)
    type = node.type
    if (typeof(type) == Interval)
        if type.type == int
            newVariable = SeaPearl.IntVar(type.start_value, type.end_value, node.id, trailer)
            SeaPearl.addVariable!(m, newVariable)
            interpreter.GLOBAL_VARIABLE[node.id] = newVariable
        elseif type.type == float
            error("Float are not permited")
        end
    elseif (typeof(type) == Domain)
        error("Domain variable are not permited")
    elseif (typeof(type) == BasicType)
        if (type.name == int)
            newVariable = SeaPearl.IntVar(0, typemax(Int64), node.id, trailer)
            SeaPearl.addVariable!(m, newVariable)
            interpreter.GLOBAL_SCOPE[node.id] = newVariable
        elseif (type.name == bool)
            newVariable =SeaPearl.BoolVar(node.id, trailer)
            SeaPearl.addVariable!(m, newVariable)
            interpreter.GLOBAL_VARIABLE[node.id] = newVariable
        else
            error("Float are not permited")
        end
    elseif (typeof(type) == ArrayVarType)
        variable_type = type.type
        start_value = type.range.start_value
        end_value = type.range.end_value
        arrayVariables = []
        for i in start_value:end_value
            if (typeof(variable_type) == Interval)
                if type.type == int
                    variable_name = node.annotations_values.values[i].value
                    push!(arrayVariables,  variable_name)
                elseif variable_type == float
                    error("Float are not permited")
                end
            elseif (typeof(variable_type) == Domain)
                error("Domain variable are not permited")
            elseif (typeof(variable_type) == BasicType)
                if (variable_type.name == int)
                    variable_name = node.annotations_values.values[i].value
                    push!(arrayVariables,  variable_name)
                elseif (variable_type.name == bool)
                    variable_name = node.annotations_values.values[i].value
                    push!(arrayVariables,  variable_name)
                else
                    error("Float are not permited")
                end
            end
        end
        interpreter.GLOBAL_VARIABLE[node.id] = arrayVariables
    end
end

function create_parameter(interpreter::Interpreter, parameter)
    if (typeof(parameter.type) == ArrayParType)
        interpreter.GLOBAL_PARAMETER[parameter.id] = parameter.expression.values
    end
end


function create_constraint(interpreter::Interpreter, constraint, trailer, m)
    if (occursin("all_different", constraint.id))
        variables = SeaPearl.AbstractIntVar[]
        #new_constraint = SeaPearl.AllDifferent(temps, trailer)
        for var in interpreter.GLOBAL_VARIABLE[constraint.expressions[1].value]
            push!(variables, interpreter.GLOBAL_VARIABLE[var])
        end
        new_constraint = SeaPearl.AllDifferent(variables, trailer)
        SeaPearl.addConstraint!(m, new_constraint)
        #push!(m.constraints, new_constraint)
        push!(interpreter.GLOBAL_CONSTRAINT, new_constraint)
    
    elseif "int_lin_eq" == constraint.id
        if (typeof(constraint.expressions[1]) == BasicExpr)
            multiplicator = interpreter.GLOBAL_PARAMETER[constraint.expressions[1].value]
        else
            multiplicator = constraint.expressions[1].values
        end
        
        variables_names = constraint.expressions[2].values
        numbers = []
        variables = []
        variablesAssignees = SeaPearl.IntVarView[]
        for i in 1:length(multiplicator)
            push!(numbers, multiplicator[i].value)
            push!(variables, interpreter.GLOBAL_VARIABLE[variables_names[i].value])
            if (numbers[i] > 0)
                push!(variablesAssignees, SeaPearl.IntVarViewMul(variables[i], numbers[i], variables_names[i].value*"_view"))
                SeaPearl.addVariable!(m, variablesAssignees[i])
                interpreter.GLOBAL_VARIABLE[variables_names[i].value*"_view"] = variablesAssignees[i]
            elseif (numbers[i] == -1)
                push!(variablesAssignees, SeaPearl.IntVarViewOpposite(variables[i], variables_names[i].value*"_view"))
                SeaPearl.addVariable!(m, variablesAssignees[i])
                interpreter.GLOBAL_VARIABLE[variables_names[i].value*"_view"] = variablesAssignees[i]
            elseif (numbers[i] < -1)
                temps = SeaPearl.IntVarViewMul(variables[i], numbers[i], variables_names[i].value*"_view")
                push!(variablesAssignees, SeaPearl.IntVarViewOpposite(temps, variables_names[i].value*"_view"))
                SeaPearl.addVariable!(m, variablesAssignees[i])
                interpreter.GLOBAL_VARIABLE[variables_names[i].value*"_view_neg"] = variablesAssignees[i]
            end
        end
        new_constraint_greater = SeaPearl.SumGreaterThan(variablesAssignees, constraint.expressions[3].value, trailer)
        new_constraint_lesser = SeaPearl.SumLessThan(variablesAssignees, constraint.expressions[3].value, trailer)
        SeaPearl.addConstraint!(m, new_constraint_greater)
        SeaPearl.addConstraint!(m, new_constraint_lesser)
        push!(interpreter.GLOBAL_CONSTRAINT, new_constraint_greater)
        push!(interpreter.GLOBAL_CONSTRAINT, new_constraint_lesser)

    elseif "int_le" == constraint.id
        #var int: a, var int: b
        #Constrains a to be less than or equal to b
    elseif "int_eq_reif" == constraint.id
        #var int: a, var int: b, var bool: r)
        #Constrains ( a = b ) ↔ r

    elseif "array_int_element" == constraint.id
        #=var int: b,
        array [int] of int: as,
        var int: c
        Constrains as [ b ] = c=#
    elseif "array_int_maximum" == constraint.id 
        #=var int: m, 
        array [int] of var int: x
        Constrains m to be the maximum value of the (non-empty) array x=#
    elseif "array_int_minimum" == constraint.id
        #=var int: m, array [int] of var int: x)
        Constrains m to be the minimum value of the (non-empty) array x=#
    
    elseif "array_var_int_element"== constraint.id
    #=var int: b,
        array [int] of var int: as,
        var int: c)
        Constrains as [ b ] = c=#
    elseif "int_abs"== constraint.id
        #var int: a, var int: b
        #Constrains b to be the absolute value of a

    elseif "int_div"== constraint.id
        #var int: a, var int: b, var int: c
        #Constrains a / b = c
    elseif "int_eq"== constraint.id
        #var int: a, var int: b)
        #Constrains a to be equal to b

    elseif "int_le_reif"== constraint.id
        #=var int: a, var int: b, var bool: r
        Constrains ( a ≤ b ) ↔ r=#
    elseif "int_lin_eq"== constraint.id
        #=array [int] of int: as,
        array [int] of var int: bs,
        int: c      
        Constrains c=∑ias[i]∗bs[i]=#
    elseif "int_lin_eq_reif"== constraint.id
        #=array [int] of int: as,
        array [int] of var int: bs,
        int: c,
        var bool: r

        Constrains r↔(c=∑ias[i]∗bs[i])=#
    elseif "int_lin_le"== constraint.id
        #=array [int] of int: as,
        array [int] of var int: bs,
        int: c
        Constrains ∑ as [ i ]* bs [ i ] ≤ c=#
    elseif "int_lin_le_reif"== constraint.id
        #=array [int] of int: as,
        array [int] of var int: bs,
        int: c,
        var bool: r
        Constrains r ↔ (∑ as [ i ]* bs [ i ] ≤ c )=#
    elseif "int_lin_ne"== constraint.id
        #=array [int] of int: as,
        array [int] of var int: bs,
        int: c
        Constrains c≠∑ias[i]∗bs[i]=#
    elseif "int_lin_ne_reif"== constraint.id
        #=array [int] of int: as,
        array [int] of var int: bs,
        int: c,
        var bool: r
        Constrains r↔(c≠∑ias[i]∗bs[i])=#
    elseif "int_lt"== constraint.id
        #var int: a, var int: b
        #Constrains a < b
    elseif "int_lt_reif"== constraint.id
        #var int: a, var int: b, var bool: r
        #Constrains r ↔ ( a < b )
    elseif "int_max"== constraint.id
        #var int: a, var int: b, var int: c
        #Constrains max( a , b ) = c
    elseif "int_min"== constraint.id
        #var int: a, var int: b, var int: c
        #Constrains min( a , b ) = c
    elseif "int_mod"== constraint.id
        #var int: a, var int: b, var int: c
        #Constrains a % b = c
    elseif "int_ne"== constraint.id
        #var int: a, var int: b
        #Constrains a ≠ b
    elseif "int_ne_reif"== constraint.id
        #var int: a, var int: b, var bool: r
        #r ↔ ( a ≠ b 
    elseif "int_plus"== constraint.id
        #var int: a, var int: b, var int: c
        #Constrains a + b = c
    elseif "int_pow"== constraint.id
        #var int: x, var int: y, var int: z
        #Constrains z = xy
    elseif "int_times"== constraint.id
        #var int: a, var int: b, var int: c
        #Constrains a * b = c
    elseif "set_in"== constraint.id
        #var int: x, set of int: S
        #Constrains x ∈ S
    elseif "array_bool_and"== constraint.id
        #array [int] of var bool: as, var bool: r
        #Constrains r↔⋀i as[i]

    elseif "array_bool_element"== constraint.id
        #var int: b,
        #array [int] of bool: as,
        #var bool: c
        #Constrains as [ b ] = c

    elseif "array_bool_or"== constraint.id
        #array [int] of var bool: as, var bool: r
        #Constrains r↔⋁ias[i]
    elseif "array_bool_xor"== constraint.id
        #array [int] of var bool: as)
        #Constrains ⊕i as[i]
    elseif "array_var_bool_element"== constraint.id
        #var int: b,
        #array [int] of var bool: as,
        #var bool: c)
        #Constrains as [ b ] = c
    elseif "bool2int"== constraint.id
        #var bool: a, var int: b
        #Constrains b∈{0,1} and a↔b=1
    elseif "bool_and"== constraint.id
        #var bool: a, var bool: b, var bool: r
        #Constrains r↔a∧b
    elseif "bool_clause"== constraint.id
        #array [int] of var bool: as,
        #array [int] of var bool: bs)
        #Constrains ⋁ias[i]∨⋁j¬bs[j]
    elseif "bool_eq"== constraint.id
        #var bool: a, var bool: b
        #Constrains a = b
    elseif "bool_eq_reif"== constraint.id
        #var bool: a, var bool: b, var bool: r)    
        #Constrains r ↔ ( a = b )
    elseif "bool_le"== constraint.id
        #var bool: a, var bool: b
        #Constrains a ≤ b
    elseif "bool_le_reif"== constraint.id
        #var bool: a, var bool: b, var bool: r)
        #Constrains r ↔ ( a ≤ b )
    elseif "bool_lin_eq"== constraint.id
        #array [int] of int: as,
        #array [int] of var bool: bs,
        #var int: c
        #Constrains c=∑ias[i]∗bs[i]
    elseif "bool_lin_le"== constraint.id
        #array [int] of int: as,
        #array [int] of var bool: bs,
        #int: c
        #Constrains ∑ias[i]∗bs[i]≤c
    elseif "bool_lt"== constraint.id
        #var bool: a, var bool: b
        #Constrains a < b
    elseif "bool_lt_reif"== constraint.id
        #var bool: a, var bool: b, var bool: r
        #Constrains r ↔ ( a < b )
    elseif "bool_not"== constraint.id
        #var bool: a, var bool: b
        #Constrains a ≠ b
    elseif "bool_or"== constraint.id
        #var bool: a, var bool: b, var bool: r
        #Constrains r↔a∨b
    elseif "bool_xor"== constraint.id
        #var bool: a, var bool: b, var bool: r
        #Constrains r↔a⊕b
    elseif "bool_xor"== constraint.id
        #var bool: a, var bool: b
        #Constrains a ⊕ b
    elseif "array_set_element"== constraint.id
        #var int: b,
        #array [int] of set of int: as,
        #var set of int: c
        #Constrains as [ b ] = c
    elseif "array_var_set_element"== constraint.id
        #var int: b,
        #array [int] of var set of int: as,
        #var set of int: c
        #Constrains as [ b ] = c
    elseif "set_card"== constraint.id
        #var set of int: S, var int: x
        #Constrains x = | S |
    elseif "set_diff"== constraint.id
        #var set of int: x,
        #var set of int: y,
        #var set of int: r
        #Constrains r = x ∖ y
    elseif "set_eq"== constraint.id
        #var set of int: x, var set of int: y
        #Constrains x = y
    elseif "set_eq_reif"== constraint.id
        #var set of int: x,
        #var set of int: y,
        #var bool: r
        #Constrains r ↔ ( x = y )
    elseif "set_in"== constraint.id
        #var int: x, var set of int: S
        #Constrains x ∈ S
    elseif "set_in_reif"== constraint.id
        #var int: x, set of int: S, var bool: r
        #Constrains r↔(x∈S)
    elseif "set_in_reif"== constraint.id
        #var int: x, var set of int: S, var bool: r
        #Constrains r↔(x∈S)
    elseif "set_intersect"== constraint.id
        #var set of int: x,
        #var set of int: y,
        #var set of int: r
        #Constrains r = x ∩ y
    elseif "set_le"== constraint.id
        #var set of int: x, var set of int: y
        #Constrains x ≤ y (lexicographic order of the sorted lists of elements)
    elseif "set_le_reif"== constraint.id
        #var set of int: x,
        #var set of int: y,
        #var bool: r
        #Constrains r↔(x≤y) (lexicographic order of the sorted lists of elements)
    elseif "set_lt"== constraint.id
        #var set of int: x, var set of int: y
        #Constrains x < y (lexicographic order of the sorted lists of elements)
    elseif "set_lt_reif"== constraint.id
        #var set of int: x,
        #var set of int: y,
        #var bool: r
        #Constrains r↔(x<y) (lexicographic order of the sorted lists of elements)
    elseif "set_ne"== constraint.id
        #var set of int: x, var set of int: y
        #Constrains x ≠ y
    elseif "set_ne_reif"== constraint.id
        #var set of int: x,
        #var set of int: y,
        #var bool: r
        #Constrains r ↔ ( x ≠ y )
    elseif "set_subset"== constraint.id
        #var set of int: x, var set of int: y
        #Constrains x ⊆ y
    elseif "set_subset_reif"== constraint.id
        #var set of int: x,
        #var set of int: y,
        #var bool: r
        #Constrains r↔(x⊆y)
    elseif "set_superset"== constraint.id
        #var set of int: x, var set of int: y)
        #Constrains x ⊇ y
    elseif "set_superset_reif"== constraint.id
        #var set of int: x,
        #var set of int: y,
        #var bool: r
        #Constrains r↔(x⊆y)
    elseif "set_symdiff"== constraint.id
        #var set of int: x,
        #var set of int: y,
        #var set of int: r
        #Constrains r to be the symmetric difference of x and y
    elseif "set_union"== constraint.id
        #var set of int: x,
        #var set of int: y,
        #var set of int: r
        #Constrains r = x ∪ y
    elseif "array_int_maximum"== constraint.id
        #var int: m, array [int] of var int: x
        #Constrains m to be the maximum value in array x .
    elseif "array_int_minimum"== constraint.id
        #var int: m, array [int] of var int: x
        #Constrains m to be the minimum value in array x .
    elseif "bool_clause_reif"== constraint.id
        #array [int] of var bool: as,
        #array [int] of var bool: bs,
        #var bool: b
        #Reified clause constraint. Constrains b↔⋁ias[i]∨⋁j¬bs[j]
    end
end


function create_solve(interpreter::Interpreter, solve_node, trailer, m)
    if (typeof(solve_node[1]) == Minimize)
        SeaPearl.addObjective!(m, interpreter.GLOBAL_VARIABLE[solve_node[1].expressions.value]) 
       #= variableToMinimize = 
        negativeVariable = SeaPearl.IntVar(-variableToMinimize.domain.max.value, -variableToMinimize.domain.min.value, "optimization",trailer)
        SeaPearl.addVariable!(m, negativeVariable)
        interpreter.GLOBAL_VARIABLE["optimization"] = negativeVariable
        SeaPearl.addObjective!(m,variableToMinimize) 
        constraint = SeaPearl.Absolute(negativeVariable,variableToMinimize, trailer)
        push!(m.constraints, constraint)
        push!(interpreter.GLOBAL_CONSTRAINT, constraint)=#
    elseif (typeof(solve_node[1]) == Maximize)
        variableToMaximize = interpreter.GLOBAL_VARIABLE[solve_node[1].expressions.value]
        negativeVariable = SeaPearl.IntVar(-variableToMaximize.domain.max.value, -variableToMaximize.domain.min.value, "optimization",trailer)
        SeaPearl.addVariable!(m, negativeVariable)
        interpreter.GLOBAL_VARIABLE["optimization"] = negativeVariable
        constraint = SeaPearl.Absolute(negativeVariable, variableToMaximize, trailer)
        SeaPearl.addObjective!(m, negativeVariable) 
        SeaPearl.addConstraint!(m, constraint)
        push!(interpreter.GLOBAL_CONSTRAINT, constraint)
        interpreter.MINIMIZATION = false
    end

end

function create_output(interpreter::Interpreter, m, node)
    output = ""
    goodSolution = findMinSolution(m.statistics.objectives)
    solution = m.statistics.solutions[goodSolution[1]]
    for variable in node
        if (variable.annotations !== nothing)
            anns = variable.annotations.annotationsList
            for ann in anns
                if (ann.id == "output_array")
                    output = output*variable.id*" = array"*string(length(ann.value[1].values))*"d("
                    for val in ann.value[1].values
                        output = output*string(val.value.value.start)*".."*string(val.value.value.stop)*", "
                    end
                    #mesvariables = array2d(1..3, 1..3, [9, 8, 7, 2, 6, 5, 4, 3, 1]);
                    output = output*"["*string(solution[interpreter.GLOBAL_VARIABLE[variable.id][1]])
                    for var in interpreter.GLOBAL_VARIABLE[variable.id][2:end]
                        output = output*", "*string(solution[var])
                    end
                    output = output*"]);\n"
                elseif ann.id == "output_var"
                    output = output*variable.id*" = " *string(solution[variable.id]) *";\n"
                end
            end
        end
    end
    output = output*"----------\n=========="
    interpreter.output = output
end

function findMinSolution(objectif)
    index = 1
    value = typemax(Int64)
    for val in 1:length(objectif)
        if objectif[val] !== nothing && objectif[val] < value
            value = objectif[val]
            index = val
        end
    end
    return (index, value)
end

function create_model(model)
    lexer = Lexer(model)
    parser = Parser(lexer)
    node = read_model(parser)
    interpreter = Interpreter(node)
    trailer = SeaPearl.Trailer()
    m = SeaPearl.CPModel(trailer)
    for variable in node.variables
        create_variable(interpreter, variable ,trailer, m)
    end
    for parameter in node.parameters
        create_parameter(interpreter,parameter)
    end
    for constraint in node.constraints
        create_constraint(interpreter, constraint, trailer, m)
    end
    
    create_solve(interpreter, node.solves, trailer, m)


    variableSelection = SeaPearl.MinDomainVariableSelection{false}()
    SeaPearl.solve!(m; variableHeuristic=variableSelection,)
    create_output(interpreter, m, node.variables)
    return interpreter
end



