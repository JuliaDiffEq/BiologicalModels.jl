### Expression Escaping ###

# Function that handles variable interpolation.
function esc_dollars!(ex)
    # If we do not have an expression: recursion has finished and we return the input.
    (ex isa Expr) || (return ex)

    # If we have encountered an interpolation, perform the appropriate modification, else recur.
    if ex.head == :$
        return esc(:($(ex.args[1])))
    else
        for i in eachindex(ex.args)
            ex.args[i] = esc_dollars!(ex.args[i])
        end
    end
    ex
end

# Checks if an expression is an escaped expression (e.g. on the form `$(Expr(:escape, :Y))`)
function is_escaped_expr(expr)
    return (expr isa Expr) && (expr.head == :escape) && (length(expr.args) == 1)
end

### Parameters/Species/Variables Symbols Correctness Checking ###

# Throws an error when a forbidden symbol is used.
function forbidden_symbol_check(syms)
    used_forbidden_syms = intersect(forbidden_symbols_error, syms)
    isempty(used_forbidden_syms) ||
        error("The following symbol(s) are used as species or parameters: $used_forbidden_syms, this is not permitted.")
end


### Catalyst-specific Expressions Manipulation ###

# Many option inputs can be on a form `@option input` or `@option begin ... end`. In both these
# cases we want to retrieve the third argument in the option expression. Further more, we wish
# to throw an error if there is more inputs (suggesting e.g. multiple inputs on a single line).
# Note that there are only some options for which we wish to make this check.
function get_block_option(expr)
    (length(expr.args) > 3) &&
        error("An option input ($expr) is missformatted. Potentially, it has multiple inputs on a single lines, and these should be split across multiple lines using a `begin ... end` block.")
    (length(expr.args) < 3) &&
        error("An option input ($expr) is missformatted. It seems that it has no inputs, which is expected.")
    return expr.args[3]
end

# Some options takes input on form that is either `@option ...` or `@option begin ... end`.
# This transforms input of the latter form to the former (with only one line in the `begin ... end` block)
function option_block_form(expr)
    (expr.head == :block) && return expr
    return Expr(:block, expr)
end

# In variable/species/parameters on the forms like:
# X
# X = 1.0
# X, [metadata=true]
# X = 1.0, [metadata=true]
# X(t)
# X(t) = 1.0
# X(t), [metadata=true]
# X(t) = 1.0, [metadata=true]
# Finds the: Variable name (X), Independent variable name(s) ([t]), default value (2.0), and metadata (:([metadata=true])).
# If a field does not exist (e.g. independent variable in `X, [metadata=true]`), gives nothing.
# The independent variables are given as a vector (empty if none given).
# Does not support e.g. "X [metadata=true]" (when metadata does not have a comma before).
function find_varinfo_in_declaration(expr)
    # Handles the $(Expr(:escape, :Y)) case:
    is_escaped_expr(expr) && (return find_varinfo_in_declaration(expr.args[1]))

    # Case: X
    (expr isa Symbol) && (return expr, [], nothing, nothing)
    # Case: X(t)
    (expr.head == :call) && (return expr.args[1], expr.args[2:end], nothing, nothing)
    if expr.head == :(=)
        # Case: X = 1.0
        (expr.args[1] isa Symbol) && (return expr.args[1], [], expr.args[2], nothing)
        # Case: X(t) = 1.0
        (expr.args[1].head == :call) &&
            (return expr.args[1].args[1], expr.args[1].args[2:end], expr.args[2].args[1],
            nothing)
    end
    if expr.head == :tuple
        # Case: X, [metadata=true]
        (expr.args[1] isa Symbol) && (return expr.args[1], [], nothing, expr.args[2])
        # Case: X(t), [metadata=true]
        (expr.args[1].head == :call) &&
            (return expr.args[1].args[1], expr.args[1].args[2:end], nothing, expr.args[2])
        if expr.args[1].head == :(=)
            # Case: X = 1.0, [metadata=true]
            (expr.args[1].args[1] isa Symbol) &&
                (return expr.args[1].args[1], [], expr.args[1].args[2], expr.args[2])
            # Case: X(t) = 1.0, [metadata=true]
            (expr.args[1].args[1].head == :call) &&
                (return expr.args[1].args[1].args[1], expr.args[1].args[1].args[2:end],
                expr.args[1].args[2].args[1], expr.args[2])
        end
    end
    error("Unable to detect the variable declared in expression: $expr.")
end

# Converts an expression of the forms:
# X
# X = 1.0
# X, [metadata=true]
# X = 1.0, [metadata=true]
# To the form:
# X(t)
# X(t) = 1.0
# X(t), [metadata=true]
# X(t) = 1.0, [metadata=true]
# (In this example the independent variable :t was inserted).
# Here, the iv is a iv_expr, which can be anything, which is inserted
function insert_independent_variable(expr_in, iv_expr)
    # If expr is a symbol, just attach the iv. If not we have to create a new expr and mutate it.
    # Because Symbols (a possible input) cannot be mutated, this function cannot mutate the input
    # (would have been easier if Expr input was guaranteed).
    (expr_in isa Symbol) && (return Expr(:call, expr_in, iv_expr))
    expr = deepcopy(expr_in)

    # Loops through possible cases.
    if expr.head == :(=)
        # Case: :(X = 1.0)
        expr.args[1] = Expr(:call, expr.args[1], iv_expr)
    elseif expr.head == :tuple
        if expr.args[1] isa Symbol
            # Case: :(X, [metadata=true])
            expr.args[1] = Expr(:call, expr.args[1], iv_expr)
        elseif (expr.args[1].head == :(=)) && (expr.args[1].args[1] isa Symbol)
            # Case: :(X = 1.0, [metadata=true])
            expr.args[1].args[1] = Expr(:call, expr.args[1].args[1], iv_expr)
        end
    end
    return expr
end
