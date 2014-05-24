r"""
All modules in this directory (and subdirectories) are devoted to a 
new implementation of inexact objects (e.g. p-adics, formal series)
in Sage.

WARNING:

    Everything is written in Python and very very slow. Do not use 
    this package for performance!


INSTALLATION::

Just copy the package in some directory /my/directory/precision
and add the three following lines to $HOME/.sage/init.sage::

    import sys
    sys.path.append("/my/directory/")
    from precision.all import *


WORKING WITH p-ADIC NUMBERS::

We first create the p-adic field Q_5::

    sage: K = QQp(5)
    sage: K
    5-adic Field

Elements in K are created by the usual syntax::

    sage: x = K(0)
    sage: x
    0
    sage: x.parent()
    5-adic Field

Note that the element `x` is exactly known (i.e. it does not carry 
a precision - or, more exactly, it carries an exact precision as we 
shall see later)::

    sage: x.is_exact()
    True

Of course, it is also possible to create inexact objects.
The simplest way to do this is to specify the precision when the
object is created::

    sage: y = K(1,20)
    sage: y
    1 + O(5^20)
    sage: y.is_exact()
    False

Another possibility to create inexact p-adics is to use the
method ``add_bigoh`` as follows ::

    sage: z = K(27); z
    27
    sage: z.add_bigoh(50)
    27 + O(5^50)

Of course, we can do basic operations with p-adics::

    sage: x+y
    1 + O(5^20)
    sage: x*y
    0
    sage: x+y-y
    O(5^20)
    sage: y^25
    1 + O(5^22)

Other usual functions are availiable::

    sage: x = K(6,20); x
    6 + O(5^20)
    sage: y = log(x); y
    9146849050361*5^1 + O(5^20)
    sage: exp(y)
    6 + O(5^20)
    sage: x.teichmuller()
    1


LAZY p-ADICS::

It may happen that computing a value requires an infinite amount
of computation. For example, if `x` is `5`-adic number `6` known
with infinite precision, computing `\log(x)` requires to expand 
the entire series defining it.

In this case, a lazy p-adics is created::

    sage: x = K(6); x
    sage: y = log(x); y
    log(6)

The term ``lazy p-adics`` means that `y` is not stored as usual
as an actual value but as a method to compute it to arbitrary 
precision. We can know whether or not a p-adic number is lazy 
by using the method ``is_evaluted``::

    sage: y.is_evaluated()
    False

By default (we will see later how to change this), a lazy p-adics
is evaluated as soon as it is possible. For instance if we ask Sage
to compute `y` to some finite precision, then the obtained result 
is an actual value::

    sage: y.add_bigoh(20)
    9146849050361*5^1 + O(5^20)

Note that `y` is nevertheless still lazy::

    sage: y
    log(6)

The same happens in the following situation (note that we have a
better precision here since `y` has valuation `1`)::

    sage: y * K(1,20)
    66367308034736*5^1 + O(5^21)

Here is a little more sophisticated example::

    sage: t = exp(log(K(126)) + 2*log(K(81)))
    sage: t   # Sage doesn't do simplifications
    exp(log(126) + 2*log(81))
    sage: t.add_bigoh(100)
    826686 + O(5^100)
    sage: 126 * 81^2
    826686

With lazy p-adics, checking equality to `0` may never terminate.
Here is a very simple example of that::

    sage: z = y-y; z
    log(6) - log(6)
    sage: z.is_zero()
    Traceback (most recent call last)
    ...
    KeyboardInterrupt

Indeed, Sage tries to compute `z` with more and more precision
since, at any finite level of precision, it cannot decide if `z`
is `0` or not.
If we want to stop the computation at some point, we can use the
parameter ``lazylimit``::

    sage: z.is_zero(lazylimit=1000)
    True

Here Sage has checked that `z` is zero up to precision `5^1000`
and has then considered that it is enough to conclude that `z` 
is indeed zero.


POLYNOMIALS AND MATRICES::

Polynomials and matrices over inexact rings (i.e. essentialy
p-adic rings and fields) are also provided by this package (but
only very basic functionalities are implemented yet).

Here is a small demo with polynomials::

    sage: P.<x> = PolynomialRing(K)
    sage: P
    Polynomial Ring in x over 5-adic Field
    sage: x
    x

    sage: y = x + K(5,20); y
    x + (5 + O(5^20))
    sage: y^2
    x^2 + (2*5^1 + O(5^20))*x + (1*5^2 + O(5^21))

    # Newton slopes
    sage: y.newton_slopes()
    [1]
    sage: (y^5).newton_slopes()
    [1, 1, 1, 1, 1]
    sage: z = x^3 + K(5,20)*x; z
    sage: z = x^2 + K(5,20); z
    x^3 + (5 + O(5^20))*x
    sage: z.newton_slopes()
    [1/2, 1/2]
    sage: (y*z).newton_slopes()
    [1/2, 1/2, 1]

    # With lazy p-adics
    sage: a = log(K(6))
    sage: y = (x + a)^3; y
    (x + log(6))^3
    sage: y.add_bigoh(5)
    (1 + O(5^5))*x^3 + (458*5^1 + O(5^5))*x^2 + (88*5^2 + O(5^5))*x + (6*5^3 + O(5^5))

Now a small demo with matrices::

    sage: MS = MatrixSpace(K,3)
    sage: MS
    MatrixSpace of 3 by 3 matrices over 5-adic Field

    sage: MS.identity_matrix()
    [ 1  0  0 ]
    [ 0  1  0 ]
    [ 0  0  1 ]

    sage: M = MS([K(i,i+5) for i in range(9)])
    sage: M
    [   O(5^5)     1 + O(5^6)     2 + O(5^7)    ]
    [ 3 + O(5^8)   4 + O(5^9)   1*5^1 + O(5^10) ]
    [ 6 + O(5^11)  7 + O(5^12)    8 + O(5^13)   ]
    sage: M*M   # look at precisions
    [ 3*5^1 + O(5^6)    18 + O(5^5)    21 + O(5^5)  ]
    [  42 + O(5^5)      54 + O(5^6)    66 + O(5^7)  ]
    [  69 + O(5^5)    18*5^1 + O(5^6)  111 + O(5^7) ]

    # With lazy p-adics
    sage: M = MS([10,2,3,log(K(6)),4,5,1,1,3])
    sage: M
    [ 2*5^1   2    3   ]
    [ log(6)  4  1*5^1 ]
    [   1     1    3   ]
    sage: N = (M*M - MS.identity_matrix())^7; N
    / [ 2*5^1   2    3   ]   [ 2*5^1   2    3   ]   [ -1  0   0  ] \^7
    | [ log(6)  4  1*5^1 ] * [ log(6)  4  1*5^1 ] + [ 0   -1  0  ] |
    \ [   1     1    3   ]   [   1     1    3   ]   [ 0   0   -1 ] /
    sage: N.add_bigoh(10)
    [  3612393 + O(5^10)    6238174 + O(5^10)    1678421 + O(5^10)   ]
    [ 50498*5^2 + O(5^10)  75689*5^2 + O(5^10)  211356*5^2 + O(5^10) ]
    [  2909822 + O(5^10)    5505946 + O(5^10)    6570434 + O(5^10)   ]


APPROXIMATION OBJECTS AND PRECISION OBJECTS::

Behind the scene, all inexact objects (p-adics, polynomials,
matrices) are represented as a couple consisting of

 - an approximation object (which is an instance of the class
   Approximation defined in the module precision.approximation)

 - a precision object (which is an instance of the class
   BigOh defined in the module precision.bigoh)

Both classes Approximation and Precision derive from Element,
which means that approximation objects and precision objects
are usual Sage elements (in particular, they have a parent).

    sage: K = QQp(5); K
    5-adic Field
    sage: x = K(3,20); x
    3 + O(5^20)

    sage: appx = x.approximation(); appx
    3
    sage: appx.parent()
    Approximation of 5-adic Ring/Field

    sage: precx = x.precision(); precx
    O(5^20)
    sage: precx.parent()
    Parent for BigOh in 5-adic Field

Even exact elements have their own precision::

    sage: K(0).precision()
    Exact

Approximations and precisions are both equipped with usual
operations::

    sage: y = K(2,17)
    sage: appy = y.approximation(); appy
    2
    sage: precy = y.precision(); precy
    O(5^17)

    sage: appx + appy
    1*5^1
    sage: appx * appy
    6

    sage: precx + precy     # addition of BigOh spaces
    O(5^17)
    sage: precx >> 2        # shifting precision
    O(5^18)
    sage: precy > precx     # check strict inclusion
    True
    sage: precy.is_exact() 
    False

If `z` is lazy, the ``z.approximation()`` can fail::

    sage: z = log(K(6)); z
    log(6)
    sage: z.approximation()
    Traceback (most recent call last)
    ...
    ApproximationError: unable to compute log to infinite precision

Nevertheless, it is possible to specify a (working) precision as
follows::

    sage: z.approximation(50)
    10678050577967888227052298704519111*5^1


The same scheme applies for polynomials and matrices as well,
except that precision objects are more involved.
Let us illustrate this with polynomials.

    sage: P.<x> = PolynomialRing(K)
    sage: a = K(1,10)*x^2 + K(5,5)*x + K(1,15); a
    (1 + O(5^10))*x^2 + (1*5^1 + O(5^5))*x + (1 + O(5^15))
    sage: a.precision()
    O(5^10)*x^2 + O(5^5)*x + O(5^15)

One can see that the precision of `a` is a vector of precisions 
(which allows to have a different precision for each individual
coefficient).

Actually, it is not completely true: by default, Sage does not
represent a precision of a polynomial as a vector of precisions
but as a Newton polygon.
The following example illustrates this::

    sage: b = K(1,10)*x^2 + K(5,20)*x + K(1,15); b
    (1 + O(5^10))*x^2 + (1*5^1 + O(5^13))*x + (1 + O(5^15))

One can see that we have lost precision on the second term. It
is due to the fact that the Newton polygon defined by the points
(2,10), (1,20) and (0,15) is just the line joining (2,10) and
(0,15). It then passes through the point (1,12.5) and it is the
reason why the precision on the `x`-coefficient is just `O(5^13)`.

Note::

    Representing precision with Newton polygons (instead of just
    vectors) is a good compromise between speed and sharpness.

We can access the underlying polygon by the following::

    sage: a.precision()._polygon
    Finite Newton polygon with 3 vertices: (0, 15), (1, 5), (2, 10)
    sage: b.precision()._polygon
    Finite Newton polygon with 2 vertices: (0, 15), (2, 10)

and also draw them *via* the method ``show()``.

It is actually possible to use other models for precisions but it
is not easy at all: we need (1) to create the adequate precision 
object, (2) to create the approximation object separately and (3)
to put them together into our polynomial. We can do it as follows::

    # Creation of the precision object
    sage: from precision.polynomial.bigoh import JaggedBigOh_polynomial
    sage: parentprec = P.precision(); parentprec
    Parent for BigOh in Polynomial Ring in x over 5-adic Field
    sage: prec = JaggedBigOh_polynomial(parentprec, [15,20,10]); prec
    O(5^10)*x^2 + O(5^20)*x + O(5^15)

    # Creation of the approximation object
    sage: parentapp = P.approximation(); parentapp
    Univariate Polynomial Ring in x over Approximation of 5-adic Ring/Field
    sage: app = parentapp([1,5,1]); app
    x^2 + 5^1*x + 1

    sage: c = P(app,prec); c
    (1 + O(5^10))*x^2 + (1*5^1 + O(5^20))*x + (1 + O(5^15))

The polynomial `c` we have just created carries a so-called ``jagged``
precision (that is a vector precision). The method ``precision_model``
gives the type of the precision carried by a polynomial::

    sage: a.precision_model()
    <class 'precision.polynomial.bigoh.NewtonBigOh_polynomial'>
    sage: b.precision_model()
    <class 'precision.polynomial.bigoh.NewtonBigOh_polynomial'>
    sage: c.precision_model()
    <class 'precision.polynomial.bigoh.JaggedBigOh_polynomial'>

When one add (or multiply) two polynomials endowed with precisions
having different precision model, conversions are made automatically::

    sage: d = a + b + c; d
    (3 + O(5^10))*x^2 + (3*5^1 + O(5^5))*x + (3 + O(5^15))
    sage: d.precision_model()
    <class 'precision.polynomial.bigoh.NewtonBigOh_polynomial'>

Conversion rules are encoded in a graph defined in the parent and
accessible as follows::

    sage: parentprec.models(graph=True)
    Digraph on 3 vertices

If we display this graph, we see the following: there are three
different vertices named respectively FlatInterval, Jagged and
Newton (and corresponding to available models) and arrows from
FlatInterval and Jagged to Newton (meaning that conversions are
made in these directions).

This graph can be modified by using the methods ``add_model``
and ``add_modelconversion``.


GROUPS OF VARIABLES::

This package allows to define groups of variables and define
a specific mode of computation for each group.

By default, all variables lie in a group named ``default``
(strictly speaking, it is not true; see next section)::

    sage: K = QQp(5)
    sage: x = K(0,20)
    sage: G = x.group(); G
    Group 'Default'

Here are the parameters of this group::

    sage: G.precision_mode()
    Individual
    sage: G.capped_relative()
    +Infinity
    sage: G.lazy_mode()
    Evaluate as soon as possible
    sage: G.remember_members()
    False

The parameter ``precision_mode`` can be:

  - ``Disable`` -- no precision is handled for any member
    of the group 

  - ``Individual`` -- each element of the group carries its 
     own precision

  - ``Collective`` -- the group itself handles the precision
    of all its members.
    Note: This mode is not yet implemented.

The ``capped_relative`` parameter means that every element of
this group will be automatically truncated up to this relative
precision. Of course, a value equal to +Infinity means that
there is no automatic truncation.

The parameter ``lazy_mode`` can be:

  - ``Disable`` -- laziness is disable for all variables in this 
    group

  - ``ASAP`` -- an expression is evaluates as soon as possible

  - ``Printing`` -- evaluations are done only before printing

  - ``Lazy`` -- evalutions are done only when necessary (e.g.
    if one ask to compute the valuation of an element)

The parameter ``remember_members`` is transparent: it specifies
whether or not the group should remember its members. It *must*
be set to ``True`` is the precision mode is ``Collective``!


It is possible to create a new group as follows (but in general
we won't never need to do this and shall rather use contexts as
explained later)::

    sage: from precision.grouping.group import GroupElements_inexact
    sage: my_group = GroupElements_inexact(
              label = 'My first group',
              precision_mode = 'individual',
              capped_relative = 50,
              lazy_mode = 'lazy',
              remember_members = True
          )
    sage: my_group
    Group 'My first group'

At first, the group is empty::

    sage: my_group.members()
    []

Now, we can create a new variable and put it in this group. We
have two solutions for this. Either, we specify the name of the
group when we create the variable::

    sage: x = K(1, group=my_group)
    sage: x
    1 + O(5^50)
    sage: x.group()
    Group 'My first group'

or, we first create the variable and the move it to the group::

    sage: y = K(5)
    sage: y = y.change_group(my_group)
    sage: y
    1*5^1 + O(5^51)
    sage: y.group()
    Group 'My first group'

In both cases, we see that the value has been truncated to
relative precision `50`.

And we can check that our group contains now two variables::

    sage: my_group.members()
    [1 + O(5^50), 1*5^1 + O(5^51)]

Note::

    When some variable is garbage collected, the group knows
    it and discards it automatically.

        sage: z = K(1, group=my_group)
        sage: len(my_group.members())
        3
        sage: z = 0
        sage: len(my_group.members())
        2

When we add or multiply (or anything else) two variables in
the same group, the result lies in the same group::

    sage: z = x + y; z
    1 + 1*5^1 + O(5^50)
    sage: z.group()
    Group 'My first group'

By the way, we can see on the above example that behaviour
of the lazy mode ``Lazy``: the value of the sum has not been 
computed! If we really want to evaluate `z`, we have to use 
the method ``evaluate()`` as follows::

   sage: z.evaluate()
   6 + O(5^50)

On the contrary, when we add two variables lying in different 
groups, an error is raised::

   sage: a = K(2)
   sage: a.group()
   Group 'Default'
   sage: a + x
   Traceback (most recent call last)
   ...
   GroupError: Can add only two elements in the same group


CONTEXT::

A context is a something like a default group but may be
local to a function or a specific part of code.

Let us first examine the global context. By default, it
is the group 'Default' we have already seen but it can
be modified as follows::

    sage: set_default_context(
              label = 'Capped Relative 20',
              precision_mode = 'individual',
              capped_relative = 20,
              lazy_mode = 'disable',
              remember_members = False
          )

Note::

    It is possible to omit some parameters. In that case,
    they are set according to the previous default context.

The previous command creates a new group labelled ``Capped
Relative 20`` and having the specified parameters and tells
Sage that, from now, any new inexact variable should belong
to this group (unless something else explicitly specified
of course)::

    sage: K = QQp(5)
    sage: x = K(1)
    sage: x    # x is automatically truncated
    1 + O(5^20)
    sage: x.group()
    Group 'Capped Relative 20'

We can go back to the original default context by using
the command ``unset_default_context()``::

    sage: unset_default_context()
    sage: x = K(1); x
    1
    sage: x.group()
    Group 'Default'


Now, let us explain what are local contexts. Imagine, we
want to write a function that has to work with specific
parameters (e.g. we want to disable locally laziness or 
precision track for speed). We can simply do it using a 
decorator::

    sage: @context(lazy_mode='disable')
    sage: def my_function(*args):
    ...       # do something
    ...       return answer

It works as follows: each time the function ``add_one`` 
is called, a new group with lazy mode disable will be
created, all inexact inputs in ``*args`` will be moved
to this group, all computations done by ``my_function``
will be performed within this group, and finally the
returned value ``answer`` will be move back to the
group defined by the default context.

Here is a concrete example::

    sage: @context(capped_relative=10)
    sage: def capped_to_10(x):
    ...       return x

    sage: x = K(25); x
    1*5^2
    sage: y = capped_to_10(x); y
    1*5^2 + O(5^12)
    sage: y.group()
    Group 'Default'

Note::

    This functionality will be particularly interesting
    when the precision mode ``Collective`` will be
    implemented.
"""
