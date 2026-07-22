
# Getting started

# ... with Julia (in REPL)

Getting started with Julia (in REPL)


!!! info 
    Link to the Pluto jl and direct start with  example...

Here, we shall **not introduce** Julia's syntax and concepts for which many tutorials are available on the web. 
Instead, we just wish to remind and highlight some simple (syntax) features that help to go easier around with JAC, and 
especially for occasional users from experiment or teaching. This reminder aims to lower the initial *threshold* for users 
that have been trained on other languages in the past. Here, we shall *pick up* some issues whose physics background 
is explained only later in other tutorial. Obviously, however, Julia is a very rich and powerful language with many 
features that go well beyond of what is (and will ever) needed for JAC.

In brief, JAC provides tools for performing atomic (structure) calculations of different kind and complexity, and for 
which further details are given in the tutorials below. To see anything from JAC, we shall first invoke the tools by:

```@example startJulia
using JenaAtomicCalculator
```

a line that will appear at the beginning of all subsequent tutorials. -- A first powerful and frequently needed feature 
refers to Julia's help pages or just the "?". By typing, for instance, atom or computation

```
? atom
```

```
search: atomic Atomic AtomicState AtomicCompass AtomicStructure @atomic @atomicswap @atomicreplace

Couldn't find atom
Perhaps you meant atomic, atan, acot, acos, htol, hton, ltoh, ntoh, Atomic, @atomic, pathof, atand, atanh, cat, 
match, catch, stat, acotd, acoth, add, ans, abs, abs2, acosd, acosh, acsc, all, all!, any, any!, asec, asin, axes, 
as, tan or Beam
No documentation found.

Binding atom does not exist.
```

we see, that `atom` itself is not a well-defined term in the JAC toolbox but that there exists a number of related 
terms, such as `Atomic`, `AtomicState` (two modules of JAC) and others. We shall not enter here the modular structure 
of the JAC toolbox but start much simpler with: 

```
? Orbital
```

```@docs    ; canonical=false
Orbital
```

which, apart from its formal meaning, is a particular data structure (`struct`) of JAC and which represents a 
relativistic orbital (function) including additional information that appears helpful in the given implementation. There 
are very many (say, more than 300) of such data struct's specified in the JAC toolbox, and thus quite obvious that
nobody will remember the details of all these definitions. Indeed, the "?" is the right and a powerful means to remind 
yourself and make use of these data structures whenever necessary. Special care has been taken that all data structures 
and functions/methods comes with a reasonable explanation (docstring) in order to work efficiently with JAC.

For instance, we might ask of what can be *added* to each other in JAC:
```
? add
```

```@docs    ; canonical=false
add
```

Apart from a short explanation, these docstring always tell the user (i) in which module the method is defined; 
(ii) which arguments it takes, including Julia's *multiple dispatch* feature as well as (iii) the type of the return 
value. All these information are typically relevant to the user, especially if some input or output does not behave as 
it should. Indeed, the complexity can grow quite rapidly, for instance, if we ask for help of what we can `generate`:

```
? generate
```

```@docs    ; canonical=false
generate
```

Well, this is quite a lot, and we shall explain some of these methods below; a similar or even larger output, you can 
generate by `? perform` as well as few other terms that are central to the implementation of JAC.

**Constructors & program control:$\quad$** Another frequent use of the (help) "?" concerns the data flow and control of 
almost all computations. In JAC, we often make use of (so-called) `Settings` that enable the user to overwrite default 
values or to *control* the computation to the extent, he or she wishes to have control. These `Settings` are context 
dependent and are different for each atomic property or process that can be computed by the JAC toolbox. They are defined 
in the various modules and need to be specified accordingly. For instance, to control the computation of transition 
probabilities for the (fine-structure) levels between given initial- and final-state configuration, one has to overwrite 
the (defaults) settings:

```
? PhotoEmission.Settings
```

```@docs    ; canonical=false
PhotoEmission.Settings
```

We shall meet these and (many) other settings quite often in the tutorials below. --- Beside of Julia's help features (?),
however, it is sometimes difficult to remember the right term or function name. In this case, it easy to make a <double-tab> 
after the dot (notation) or to make use of the (Unix/Linux) `grep` command within the `JenaAtomicCalculator/src` directly. Similar 
line-search commands will exist also at other platforms. In particular, for those of you who wishes to support and extend 
the JAC toolbox, the dot expansion and the `grep` command will be found very helpful, perhaps more than other search tools.


!!! info
    Users can also refer to the API Reference section that provides a comprehensive list of JAC declared `Type`s and `Function`s 
    for selected atomic processes.

**Use of constructors:$\quad$** Another Julia feature, that is frequently applied in JAC, is the successive definition of 
constructors in order to set-up complex data structures. This features is applied, for instance, in order to define an 
`Atomic.Computation` or a `Cascade.Computation` as a whole. We shall explain these rather complex data types below in 
different tutorials. The same issue appears however already at a much simpler level. For example, if we wish to select 
(specify) a number of levels from a multiplet prior to some particular -- configuration interaction -- computation, 
we can make use of a
```
? LevelSelection
```

```@docs ; canonical=false
LevelSelection
```
Apart from the logical flag `active`, such a level selection requires to either specify a list of level numbers (indices) 
or *level symmetries*

```@example startJulia
LevelSelection(true, indices= [i for i in 1:11])
```

```@example startJulia
LevelSelection(true, symmetries= [LevelSymmetry(1//2, Basics.plus)])
```

Here, we made use of a LevelSymmetry to specify the overall rotational $J^P$ symmetry of atomic levels.

```
? LevelSymmetry
```

```@docs ; canonical=false
LevelSymmetry
```

As seen from this definition, the level symmetry just comprises the total angular momentum (of type `AngularJ64`) and 
the parity of the level (of type `Parity`). Therefore, the specification of a list of level symmetries in `LevelSelection` 
already requires to nest four constructors in order make the level selection explicit: (i) For the angular momentum, 
(ii) the parity, (iii) the level symmetry and (iv) to create a list (array) of such level symmetries. All the constructors 
can be specified and built together also in subsequent steps, such as:

```@example startJulia
J1    = AngularJ64(1//2);           J2 = AngularJ64(5//2)  
pl    = Basics.plus;                mn = Basics.minus
lsym1 = LevelSymmetry(J1, pl);      lsym2 = LevelSymmetry(J2, mn)
levelsyms = [lsym1, lsym2]
```

or simply by *nesting* all the information within a single step
```@example startJulia
levelsyms = [LevelSymmetry(AngularJ64(1//2), Basics.plus), LevelSymmetry(AngularJ64(5//2), Basics.minus)]
```

Both way have their pros and cons, and often some *mixture* is applied where complex constructors are first assigned to 
some variables, and which are later utilized to built up constructors of higher complexity. --- To finally specify aǹ 
instance of a `LevelSelection`, we use (onc more) its second constructor above:

```@example startJulia
LevelSelection(true, indices=[1,2,3], symmetries=levelsyms)
```

and which will tell the JAC program to compute the lowest three levels (1, 2, 3) as well as all levels with 1/2+ and 5/2- 
symmetry. Apart from the selection of individual levels, it is often helpful for the computation of atomic processes to 
make a prior `LineSelection` and in some cases even a `PathwaySelection` as, for instance, for dielectronic recombination 
processes. Some of these features will be explained below in subsequent tutorials of JAC:

**Functions & methods:$\quad$** Like most other languages, Julia is based on the successive work through functions and 
methods; a **function** is first of all specified by its name and it maps a tuple of argument values upon a return value. 
For instance, the function

```@example startJulia
function addSomething(a, b)
    c = a + b
end
addSomething(3.1, 5//2)
```

can be used to *add* two numbers, rather independent of their particular type, and which are *infered* here automatically. 
However, additional type declarations might help to *specialize* a function and to ensure **type stability**:

```@example startJulia
function addSomething(a::Int64, b::Int64, c::Int64)
    d = a + b + c
end
```

While the function name is the same in both of these examples above, Julia carefully distinguishes between these two 
**methods** of the function `addSomething` that may differ by the type *and/or* the number of arguments. This multiple 
*use* **(dispatch)** of function name enables the user to write highly specialized code. Although a proper (and specialized) 
definition of functions is often very important for the performance of the program, we shall not discuss such technical 
issues here. Let us just mention, that a function/method may also return `nothing`:

```@example startJulia
typeof(nothing)
```

In JAC, the value `nothing` is usually returned by all *display* functions that print some selected data or tabulation 
to screen or elsewhere but does not return a value otherwise.

**Code failures:$\quad$** Beside of its large flexibility and user-friendliness, JAC might terminate from time to time 
for *non-obvious* reasons. Since JAC is first of all a *physics code*, no attempt has been made that all possible errors 
are fully captured and recovered by the program. Wrong input parameters or an inappropriate use of contructors will often 
lead to errors that cannot be resolved by the program. While some of this input can be readily recognized as wrong, and 
then lead to a proper error message, other wrong data may appear dynamically and cannot be captured with a reasonable 
overhead of the code.  In JAC, therefore, several conditional `if ... elseif ... else ... end` blocks include an additional 
clause `error("stop a")` or similar; these are clauses, which due to a first design of a function should never be entered, 
but this appears not to be true in all cases. The use of these (fully) *non-instructive* error message have still a great 
advantage due to Julia: If not switched-off explicitly, Julia always reports for all program failures the hierarchy of 
call's, that are made before the error occurs, and lists these calls together with the file and line number of source code.
For this reason, an `error("stop a")` readily shows the position where something unexpected occurs. A short inspection of 
the corresponding (line of the) source code often help to understand of what went wrong internally.

**Julia macros:$\quad$**  What can one do, if the (source) code itself does not tell so much about the problem ? --- 
In this case, it is often useful to include some additional **printouts**  near to the line in question into the code 
and to re-run it again. There are different ways (`@show`, `print()`, `println()`) to place printout in the code; 
cf. https://julialang.org/learning/  A particular quick and useful way makes use of the Julia macro:

```@example startJulia
@show levelsyms
```

which simply repeats the *names* of the variables together with their values. Of course, the values of several such 
variables can be shown within the same call:

```@example startJulia
wa = 5;   wb = [2.0, pi];   wc = ones(3)
@show wa, wb, wc
```

Indeed, this `@show` macro makes printout very easy. There are many macros (all starting with `@`) in Julia which need 
not to be considered here. We just mention that `@time` in front of a Julia command (block) will take and display the CPU 
time that is necessary to run this line(s):

```@example startJulia
@time rand(50000)   ;
```


# ... with JAC (in REPL)

Getting started with JAC (in REPL)

!!! info
    JAC user guide pdf .... link

```@example startJAC
using JenaAtomicCalculator
```

### Welcome to **JAC**, the **Jena Atomic Calculator**

... that provides various tools for performing atomic (structure) calculations of different kinds and complexities. 
Apart from the computation of atomic (many-electron) amplitudes, properties and processes, **JAC supports interactive, 
restricted-active space (RAS) and cascade computations**. It also help perform a few simple *hydrogenic* and 
*semi-empirical* estimates as well as simplify symbolic expressions from Racah's algebra. --- 
Let's first use  `? JenaAtomicCalculator` in order to obtain more information about this toolbox:

```
? JenaAtomicCalculator
```

```@docs ; canonical=false
JenaAtomicCalculator
```
H'm, this tells us a lot of details which we still need to better understand.

In the design of JAC, we first of all **aim for a precise language** that (i) is simple enough for both, seldom 
and a more frequent use of this package, (ii) highlights the underlying physics and (iii) avoids most technical 
slang that is often unnecessary but quite common to many other codes. An intuitive picture about the level or hyperfine 
structure of an atom, its properties as well as possible excitation and/or decay processes should (always) come first 
in order to generate the desired data: By making use of suitable data types (`struct`), **we indeed wish to 
introduce a language close to the underlying formalism.** --- While JAC is overall based on a rather large 
number $(> 300)$ of such types, a few simple examples are:

  + (atomic) `Shell`:                 $\quad$1s, 2s, 2p, ...
  + `Subshell`:                       $\quad$1s_1/2, 2s_1/2, 2p_1/2, 2p_3/2, ...
  + (electron) `Configuration`:       $\quad$1s^2 2s^2 2p^6 3s $\quad$  or $\quad$  [Ne] 3s, ...
  + `Level`:                          $\quad$1s^2 2s^2  ^1S_0, ...
  
and many other terms (types) that we shall explain later.  


Let us simply start, for instance, with specifying and assigning the $1s$ and $2p$ shells:

```@example startJAC
w1s = Shell("1s")
w2p = Shell("2p")
```

Similarly, we can readily specify and assign any (relativistic) subshell:

```@example startJAC
Subshell("2p_1/2"),   Subshell("2p_3/2")  
```

In JAC, we make use of these `Shell`'s and `Subshell`'s whenever they will naturally occur in describing the level 
structure or the excitation, decay or occupation of an atom, and this both at input and output. If you have 
*forgotten* how to specify such a subshell (constructor), simply *ask*:

```
? Subshell
```

```@docs ; canonical=false
Subshell
```

Of course, we can interactively also specify any electron configuration:

```
? Configuration
```

```@example startJAC
wc1 = Configuration("1s^2 2s^2 2p^5")
wc2 = Configuration("[Ar] 4s^2 3d^5")
```
!!! info 
    For specific processes users can find the list of types and functions in the API Reference

This input just shows three (very) simple examples and how the details of some computation can be readily specified 
in line with our basic understanding of the atomic shell model.

Julia (and JAC) is a very suitable and powerful framework for running -- many-electron -- atomic computations. 

Of course, there are many other features that make Julia & JAC as powerful as it is: For example, the user may pre-define 
and overwrite the **units** in which he wishes to communicate with JAC. These units determine how (most of) the input 
data are interpreted as well as output data are displayed in tabulations or to screen. The current defaults settings 
for the units can be seen by typing:

```@example startJAC
Basics.display("settings")
```

which show that energies are taken/printed in eV, rates in 1/s, etc. Apart from modifying these defaults directly in the 
source code, the can be *overwritten* by the user at any time of the program executation. This is done by means of 
the function

```
? Defaults.setDefaults
```

```@docs
Defaults.setDefaults
```
which enables one to re-define various **global values** of JAC. If we wish to enter/display energies in **Kaysers** or 
cross sections in atomic units, we can simply type:

```@example startJAC
Defaults.setDefaults("unit: energy", "Kayser")
Defaults.setDefaults("unit: cross section", "a.u.")
```

Here, again `nothing` is returned but the corresponding global constants are now changed.

```@example startJAC   
Basics.display("settings")
```
Apart from the default units, one can similarly *overwrite* the method that is use for the generation and normalization 
of continuum orbitals and several others. Although called *global*, the corresponding values can be accesses just in 
two ways. (i) The **global constants**, such as the electron mass, the speed of light, the fine-structure constant $\alpha$, 
etc., are accessed via the function:

```
? Defaults.getDefaults
```

```@docs
Defaults.getDefaults
```

```@example startJAC
Defaults.getDefaults("alpha")
Defaults.getDefaults("electron rest energy")
Defaults.getDefaults("unit: energy")
```

(ii) These **global values** are frequently applied in order to -- internally or externally -- convert physical numbers 
into units of the same dimension. This is done by the function:

```
? Defaults.convertUnits
```

```@docs
Defaults.convertUnits
```

This function is called at many places within JAC to generate tables where all physical data are printed out in the 
pre-specified units:

```@example startJAC
Defaults.convertUnits("energy: from atomic", 1.0)
```

With the given user-selection, this is equivalent to:

```@example startJAC
Defaults.convertUnits("energy: from atomic to Kayser", 1.0)
```

In JAC, the call of this function is often combined with some proper formatting of the results, such as:

```@example startJAC
using Printf
@sprintf("%.4e", Defaults.convertUnits("energy: from atomic", 1.0))
```

