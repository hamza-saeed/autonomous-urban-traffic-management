
% This file contains the domain model in Prolog form

domain_name(urbantraffic).

types([junction,link,stage,section]).

% declarations of predicates - <type>(args) means the arguments are of type <type>
predicates([
% is true if stage P is in intergreen 
     inter(stage(P)), 
% is true if stage P is on green
     active(stage(P)), 
% asserts that stage P1 follows stage P (stage is a 2 place predicate)
     next(stage(P,P1)), 
% asserts that stage P is a stage in junction I
     contains(junction(I),stage(P)),
% this is a trigger used in the events below to start intergreen in junction I
     trigger(junction(I)), 
% sections are used to regulate flow through a link - this is a pointer to section I
     current(section(S),link(L)),
% this asserts that section S1 follows S (current is a 2 place predicate)
     nextsection(section(S,S1))
     ]).

% declarations of functions which hold numerical values
functions([
% this holds the value of the flow rate during X between R1 and R2 (static)
     f(turnrate,[stage(X),link(R1),link(R2)]),
% this holds the length of the intergreen time after stage P (static)
     f(interlimit,[stage(P)]),
% this holds the number of PCUs in a link (dynamic)
     f(occupancy,[link(R)]),
% this holds the maximum number of PCUs in a link (static)
     f(capacity,[link(R)]),
% this holds the maximum length of the green time for stage P
     f(maxgreentime,[stage(P)]),
% this holds the minimum length of the green time for stage P
     f(mingreentime,[stage(P)]),
% this records the amount of time a junction (hence a stage) is on green (dynamic)
     f(greentime,[junction(I)]),
% this records the amount of time a junction (hence a stage) is in intergreen (dynamic)
     f(intertime,[junction(I)]),
% this holds the number of PCUs in a section S (dynamic)
     f(store,[section(S)]),
% this holds the number of PCUs at the end of link R (dynamic)
     f(queue,[link(R)]),
% this records the amount of time traffic in a section have been waiting in link I (dynamic)
     f(sectiontime,[link(I)]),
% this records the length of time traffic stay in a section  in link I (dynamic)
     f(maxsectiontime,[link(I)])
     ]).

event(maxgreenreached,
% parameters
     [stage(P),junction(I)],
% preconditions
     [active(P),contains(I,P),compare(>=,greentime(I),maxgreentime(P))],
% effect - predicates made true
     [trigger(I)],
% effect - predicates made false
     [],
% effect - function values changed 
     []
).

process(keepgreen,
% parameters
     [stage(P),junction(I)],
% preconditions
     [active(P),contains(I,P),compare(<,greentime(I),maxgreentime(P))],
% effects
     [increase(greentime(I),op(*,t,1))]
).

process(flowrun_green1,
% parameters
     [stage(P),link(R1,R2),section(C,C1)],
% preconditions
     [active(P),
      compare(>,queue(R1),0.0),
      compare(<,capacity(R1),9000.0),
      compare(<,capacity(R2),9000.0),
      current(C,R2),
      nextsection(C1,C),
      compare(>,turnrate(P,R1,R2),0.0),
      compare(<,occupancy(R2),op(*,0.9,capacity(R2)))],
% effects
     [increase(occupancy(R2),op(*,t,turnrate(P,R1,R2))),
      decrease(occupancy(R1),op(*,t,turnrate(P,R1,R2))),
      increase(store(C1),op(*,t,turnrate(P,R1,R2))),
      decrease(queue(R1),op(*,t,turnrate(P,R1,R2)))]
).

process(flowrun_green2,
% parameters
     [stage(P),link(R1,R2),section(C,C1)],
% preconditions
     [active(P),
      compare(>,occupancy(R1),0.0),
      compare(>,capacity(R1),9000.0),
      compare(<,capacity(R2),9000.0),
      current(C,R2),
      nextsection(C1,C),
      compare(>,turnrate(P,R1,R2),0.0),
      compare(<,occupancy(R2),op(*,0.9,capacity(R2)))],
% effects
     [increase(occupancy(R2),op(*,t,turnrate(P,R1,R2))),
      decrease(occupancy(R1),op(*,t,turnrate(P,R1,R2))),
      increase(store(C1),op(*,t,turnrate(P,R1,R2)))]
).

process(flowrun_green3,
% parameters
     [stage(P),link(R1,R2)],
% preconditions
     [active(P),
      compare(>,queue(R1),0.0),
      compare(<,capacity(R1),9000.0),
      compare(>,capacity(R2),9000.0),
      compare(>,turnrate(P,R1,R2),0.0)],
% effects
     [increase(occupancy(R2),op(*,t,turnrate(P,R1,R2))),
      decrease(occupancy(R1),op(*,t,turnrate(P,R1,R2))),
      decrease(queue(R1),op(*,t,turnrate(P,R1,R2)))]
).

process(flowrun_green4,
% parameters
     [stage(P),link(R1,R2)],
% preconditions
     [active(P),
      compare(>,occupancy(R1),0.0),
      compare(>,capacity(R1),9000.0),
      compare(>,capacity(R2),9000.0),
      compare(>,turnrate(P,R1,R2),0.0)],
% effects
     [increase(occupancy(R2),op(*,t,turnrate(P,R1,R2))),
      decrease(occupancy(R1),op(*,t,turnrate(P,R1,R2)))]
).

action(switchphase,
% parameters
     [stage(P),junction(I)],
% preconditions
     [active(P),contains(I,P),compare(>,greentime(I),mingreentime(P))],
% effect - predicates made true
     [trigger(I)],
% effect - predicates made false
     [],
% effect - function values changed 
     []
).

process(section-tick,
% parameters
     [link(R)],
% preconditions
     [compare(<,sectiontime(R),maxsectiontime(R))],
% effects
     [increase(sectiontime(R),op(*,t,1))]
).

event(move-flow-into-queue,
% parameters
     [link(R),section(C,C1)],
% preconditions
     [current(C,R),
      nextsection(C,C1),
      compare(>=,sectiontime(R),maxsectiontime(R))],
% effect - predicates made true
     [current(C1,R)],
% effect - predicates made false
     [current(C,R)],
% effect - function values changed 
     [assign(sectiontime(R),0),
      assign(queue(R),op(+,queue(R),store(C))),
      assign(store(C),0)]).

event(trigger-inter,
% parameters
     [stage(P),junction(I)],
% preconditions
     [trigger(I),
      active(P),
      contains(I,P)],
% effect - predicates made true
     [inter(P)],
% effect - predicates made false
     [trigger(I),active(P)],
% effect - function values changed 
     [assign(greentime(I),0)]).

process(keepinter,
% parameters
     [stage(P),junction(I)],
% preconditions
     [inter(P),
      contains(I,P),
      compare(<,intertime(I),interlimit(P))],
% effects
     [increase(intertime(I),op(*,t,1))]).

event(trigger-change,
% parameters
     [stage(P,P1),junction(I)],
% preconditions
     [inter(P),
      contains(I,P),
      next(P,P1),
      compare(>=,intertime(I),op(-,interlimit(P),0.1))],
% effect - predicates made true
     [active(P1)],
% effect - predicates made false
     [inter(P)],
% effect - function values changed 
     [assign(intertime(I),0)]).
