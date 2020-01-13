
%%%%%%%%%%%%% COURSEWORK 2019 AAIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% PDDL+ SIMULATOR - L McCluskey 27/11/19 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% VERSION 1.0.0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% SWI Prolog %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%  THIS VERSION HAS A BIT OF AN ANSWER...

% PRELIMINARY FACTS

% contains pre-compiled definitions of list functions ..
:- use_module(library(lists)).

% windowns or mac
computer(windows).
% windows.
% linux

% trace ON or OFF
sim_trace(on).


% WHERE INPUT COMES FROM
initial('../Users/lee/Desktop/coursesim/initial_state.pl') :- computer(mac).
% computer(windows).
initial('C:/Users/Hamza/Desktop/AAIS Assignment/AAIS/scenario5.pl') :- computer(windows).
% Domain Model File
domainmodel('../Users/lee/Desktop/coursesim/domain_model.pl') :- computer(mac).
domainmodel('C:/Users/Hamza/Desktop/AAIS Assignment/AAIS/domain_model.pl') :- computer(windows).

% WHERE to put log information
tracefile('../Users/lee/Desktop/coursesim/trace.txt') :- computer(mac).
tracefile('C:/Users/Hamza/Desktop/AAIS Assignment/AAIS/trace.txt') :-computer(windows).

%%%%% ************* INITIALISE + go   **************************************
%% Calling predicate .... Time = Simulation time, Inc = Time increment (stick to 1)
sim(Time,Inc) :-
% read in the initial state file and domain model
    initial(I),
    domainmodel(D),
    compile(I),
    compile(D),
% do some consistency checking of links
    check_links,
    check_links_max,
% call simulation
    simulate(Time,Inc).

%%%%% *******************  SIMULATION CODE ********************************

% e.g. simulate(900,1).
% Stops after Sim_Time no of clicks

% get ready for simulation loop
simulate(Sim_Time,Delta) :-
% assert initial node ready for simulation loop:
        init(_, _, _, Init_DYN, STATICS, Init_PLAN, Init_PRESSES, _),
        assert(statics(STATICS)),
        assert(node(1,0,Init_DYN,Init_PLAN,Init_PRESSES)),
% for logging status of roads:
        assert_saturation_levels,
        click_on(Sim_Time,Delta),
        !.

init:-
	initial(I),
	compile(I),
	init(_, _, _, DYN, STATICS, PLAN, _,_),
	assert(dyn(DYN)),
	assert(statics(STATICS)),
	assert(plan(PLAN)),

        %Flags to show whether inlinks have already been changed
        AlbertSouthFlag is 0,
        CallansSouthFlag is 0,
        HeptonSouthFlag is 0,
        TradesNorthFlag is 0,

        %Store flags so they are accessible anywhere in the code
        nb_setval(albertSouthFlag,AlbertSouthFlag),
        nb_setval(callansSouthFlag,CallansSouthFlag),
        nb_setval(heptonSouthFlag,HeptonSouthFlag),
        nb_setval(tradesNorthFlag,TradesNorthFlag)
        .


maximize(X,V):-
	statics(STATICS),
	member(equals(maxgreentime(X),V), STATICS),
	retract(plan(PLAN)),
	delete(PLAN,equals(defaultgreentime(X),_),REDUCED_PLAN),
	assert(plan([equals(defaultgreentime(X),V) | REDUCED_PLAN])).
read_in_state :-
		initial(I), compile(I).

% 1. test with occ_link(hepton_south, V).
occ_link(L, V) :-
		init(_, _, _, DYN, _, _, _, _),
		member( equals(occupancy(L),V),DYN).

% 2. test with junction_name_end(hepton_south, J).
junction_name_end(L, J) :-
		init(_, _, _, _, STATICS, _, _, _),
		member( equals(turnrate(Stage,L,_),_), STATICS),
		member( contains(J, Stage), STATICS).

% 3. test with stage_len(fox_stage0, SL).
stage_len(S, SL) :-
		init(_, _, _, _, _, PLAN, _, _),
		member(equals(defaultgreentime(S),SL), PLAN).

% 4. test with list_of_stages(fox, SList).
list_of_stages(J, SList) :-
		init(_, _, _, _, STATICS, _, _, _),
                get_stages(J,STATICS,  SList).
get_stages(_, [ ], [ ]).
get_stages(J, [contains(J,S)  | RestSTATICS],   [S | RestS] ) :-
		get_stages(J, RestSTATICS,  RestS),!.
get_stages(J, [ _ | RestSTATICS],   RestS ) :-
		get_stages(J, RestSTATICS,  RestS),!.

% 5. test with cyc_length(fox, CL).
cyc_length(J, CL) :-
		init(_, _, _, _, STATICS, PLAN, _, _),
		get_cyc_length(J, PLAN, STATICS, 0,    CL).
get_cyc_length(_, [], _, CL,    CL).
get_cyc_length(J, [equals(defaultgreentime(S),SL) | RS], STATICS, CLx, CL) :-
		member( contains(J, S), STATICS),
		member(equals(interlimit(S),V), STATICS),
                CLx1 is CLx + SL + V,
		get_cyc_length(J,RS,STATICS, CLx1, CL),!.
get_cyc_length(J, [_ | RestS], STATICS, CLx, CL) :-
		get_cyc_length(J, RestS, STATICS, CLx, CL),!.


% 6.  test with
% total_flow_rate(hepton_south,fox_stage1, X).
% total_flow_rate(station_west,vocation_stage0, X).

total_flow_rate(L,S, FL) :-
		init(_, _, _, _, STATICS, _, _, _),
		get_total_flow_rate(L,S,STATICS,0,  FL).

get_total_flow_rate(_,_, [], FL,    FL).
get_total_flow_rate(L,S, [ equals(turnrate(S,L,_),V) | RestS], FLx, FL) :-
                FLx1 is FLx + V,
		get_total_flow_rate(L,S, RestS, FLx1, FL),!.
get_total_flow_rate(L,S, [_ | RestS], FLx, FL) :-
		get_total_flow_rate(L,S, RestS, FLx, FL),!.


% 7. test with
% total_flow_out_cycle(station_west, FL).
% total_flow_out_cycle(hepton_south, FL).

total_flow_out_cycle(L, FL) :-
                junction_name_end(L, J),
		list_of_stages(J, SList),
                acculmulate_flow(L,SList, FL).

acculmulate_flow(L, [S | SR], FL) :-
		total_flow_rate(L,S, TFR),
                stage_len(S, STL),
                Stage_Flow is STL*TFR,
                acculmulate_flow(L,SR, FL1),
                FL is FL1 + Stage_Flow,
		!.
acculmulate_flow(_,[], 0).


% 8. test with average_flow_out(hepton_south, X).
average_flow_out(L, AV) :-
		junction_name_end(L, J),
		cyc_length(J, CL),
		total_flow_out_cycle(L, FL),
		AV is FL/CL,
                !.

% 9. Week 10 / Q3: ADJUST_PLAN
% test with:
%
adjust_default_plan :-
                init(_, _, _, _, STATICS, PLAN, _, _),
                adjust_plan(station_west, PLAN,_,STATICS,  PLAN1),
                nl,write(PLAN1).

adjust_plan(BL, PLAN,_,STATICS,  PLAN1) :-
                occ_link(BL, OV),
% so we can test it in the initial state:
                OV > -20,
                junction_name_end(BL,J),
                list_of_stages(J,SList),
                member(S,SList),
% find any stage with Flow Out > 0 - in coursework there is only one with flow out > 0
                total_flow_rate(BL,S, FL),
                FL > 0,
% rest is taken from the "maximise" solution - Week 10
                member(equals(maxgreentime(S),MV),  STATICS),
                delete(PLAN,    equals(defaultgreentime(S),_),  REDUCED_PLAN),
                PLAN1 = [equals(defaultgreentime(S),MV) | REDUCED_PLAN ] .

% if code above fails .. return same plan.
adjust_default_plan(_, PLAN,_,  PLAN).

%Method to change stage lengths
change_stage_lengths(Link, NewS1Length, NewS0Length, Delta, PLAN3) :-

    %Get the Junction Name
    junction_name_end(Link,Junction),
    %Get individual stages
    list_of_stages(Junction, StageList),
    nth0(0,StageList,Stage0),
    nth0(1,StageList,Stage1),
    %Get turnrate for Stage 1
    total_flow_rate(Link, Stage1, FlowRate),
    %Get cycle length for junction
    cyc_length(Junction, CycleLength),
    %Get flow from outside
    get_flow_from_outside(Link,FlowOut),
    %Calculate the new stage length
    X is CycleLength * FlowOut,
    NewS1Length is X / FlowRate,
    %Now calculate new stage 0 length keeping time fixed
    stage_len(Stage0,Greentime1),
    stage_len(Stage1,Greentime2),
    CycleGreentime is Greentime1 + Greentime2,
    %Absolute value
    NewS0Length is abs(CycleGreentime - NewS1Length),
    %Add/Remove delta
    NewS1LengthPlusDelta is NewS1Length + Delta,
    NewS0LengthMinusDelta is NewS0Length - Delta,
    %Call method to change stage lengths
    change_default_greentime(Stage1,NewS1LengthPlusDelta,Stage0,NewS0LengthMinusDelta,PLAN3),
    nl,write('Changed stage lengths for '),write(Link),nl,
    !.

%Returns flow from outside
get_flow_from_outside(Link,FlowFromOut) :-
   % init(_, _, _, _, STATICS, _, _, _),
    get_statics(STATICS),
    member( equals(turnrate(fake,outside,Link),FlowFromOut), STATICS).

%Changes greentime for stage 0 and stage 1
change_default_greentime(Stage1,S1NewValue, Stage0, S0NewValue,PLAN3):-
    retract(plan(PLAN)),
    %Delete previous values
    delete(PLAN,equals(defaultgreentime(Stage1),_), REDUCED_PLAN),
    delete(REDUCED_PLAN,equals(defaultgreentime(Stage0),_), REDUCED_PLAN_2),
    %Add new values
    NEW_PLAN = ([equals(defaultgreentime(Stage1),S1NewValue) | REDUCED_PLAN_2]),
    PLAN3 = ([equals(defaultgreentime(Stage0),S0NewValue) | NEW_PLAN]),
    assert(plan(PLAN3)).

%%%%% END OF SIMULATION - WRITE OUT INFORMATION / RESULTS TO SCREEN AND FILE

click_on(Sim_Time,_) :-
        node(_,Actual_Time,DYN,PLAN,_),
        Sim_Time =< Actual_Time,
        get_statics(STAT),
        write_end_stuff(STAT,PLAN,DYN,Sim_Time),
        tracefile(Trace),
        tell(Trace),nl,nl,
        write_end_stuff(STAT,PLAN,DYN,Sim_Time),
        write_out_turns(STAT),
        told,
        !.

% MAIN SIMULATION LOOP  ***********

click_on(Sim_Time,Delta) :-
%
        do_actions,
        do_events,
        do_processes(Delta),

%
% now move on time T to T1 and give node N a new number N1
        retract(node(N,T,DYN,PL,PR)),
        N1 is N+1,
        T1 is T+Delta,
        assert(node(N1,T1,DYN,PL,PR)),
        log_trace(T,N,DYN,PL,PR),

% iterate
        click_on(Sim_Time,Delta),
        !.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/* Simulaton of actions */
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_actions :-
        retract(node(N,T,DYN,PLAN,PRESS)),
        get_statics(STATICS),

%
%       if any default time is reached at J, put a 'tigger(J)' in dynamics
%       to signal that the stage must change (going through intergreen first)
        do_plan_defaults(PLAN,DYN,STATICS, Triggers),
        append(Triggers,DYN, DYN1),
%
%       if any junctions (pelicans) have been waiting till now (T) to go on, put a 'tigger(J)' in dynamics
%       implement any outside - sourced flow rate changes that are time depedent waiting to happen
        do_plan_waiting(PLAN,T, PLAN1,MoreTriggers),
        append(MoreTriggers,DYN1, DYN2),
%
%       if any pelican is pushed at this instant, record a 'waiting' change_state as
%       long as its not already on red OR it has not already been pushed and is
%       waiting to turn on
        do_presses(PRESS,PLAN1,DYN2,T,  BIT,PRESS2),
        append(PLAN1,BIT, PLAN2),

        %Monitor queues on albert south - change if over 10 PCUs

        queue_monitor(albert_south,albertSouthFlag,DYN2,PLAN2,PLAN3),

        %Monitor queues on callans south - change if over 10 PCUs
        queue_monitor(callans_south,callansSouthFlag,DYN2,PLAN3,PLAN4),

        %Monitor queues on hepton south - change if over 10 PCUs
        queue_monitor(hepton_south,heptonSouthFlag,DYN2,PLAN4,PLAN5),

        %Monitor queues on trades north - change if over 10 PCUs
        queue_monitor(trades_north,tradesNorthFlag,DYN2,PLAN5,PLAN6),

         N1 is N+1,
         assert(node(N1,T,DYN2,PLAN6,PRESS2)),

         !.

%   if any pelican is pushed at this instant, record a 'waiting' change_state
%     assume presses are supplied in ASCENDING order for efficiency
%     assume same buttun CANNOT be pressed > 1 in the same instant

do_presses([press(T,J)|Rest],PLAN,DYN,T, [change_stage(T1,J)],Rest) :-
% crossing not already activated
        \+ member(inter(J),DYN),
% not already been pushed and waiting to go on
        \+ member(change_stage(_,J),PLAN),
% not already about to change state
        \+ member(trigger(J),DYN),
% get the wait time, which must be greater than 0
        member(wait(J,WT), PLAN),
        T1 is T + WT,
        !.
% presses that are redundant
do_presses([press(T,_)|Rest],_,_,T, [],Rest) :- !.
% no need to look through all list as presses are in ascending order
do_presses(P,_,_,_,[],P) :- !.


% if (active ?p) (contains ?i ?p) (> (greentime ?i) default ?p )
% add (trigger ?i)
do_plan_defaults([equals(defaultgreentime(Stage),DT)|Rest],DYN,ST, [trigger(J)|RestT] ) :-
        member(contains(J,Stage),ST),
        member(active(Stage), DYN),
        member(equals(greentime(J),GT),DYN),
        GT >= DT,
        do_plan_defaults(Rest,DYN,ST, RestT),
        !.
do_plan_defaults([_|Rest],DYN,S, Trig) :-
        do_plan_defaults(Rest,DYN,S, Trig),
        !.
do_plan_defaults([],_,_,[]) :- !.

% iterate through PLAN and execute anything timed to happen now (that is T)
% last argument is the 'list of triggers' to add to DYN to change stages
do_plan_waiting([change_stage(T,J)|Rest],T, Rest2,[trigger(J)|RestT]) :-
        do_plan_waiting(Rest,T, Rest2,RestT),
        !.
do_plan_waiting([change_flow(T,Link_in,NewFlow)|Rest],T, Rest2,Trig) :-
% DYN, PLAN and PRESSES now all stored in NODE so no need to keep them in init ..
        member(equals(turnrate(fake,outside,Link_in),INrate),STAT),
        retract( statics(STAT) ),
        delete(STAT,equals(turnrate(fake,outside,Link_in),INrate),NEW_STAT),
        assert( statics( [equals(turnrate(fake,outside,Link_in),NewFlow) | NEW_STAT] ) ),
        tell(user),nl,write(" CHANGED FLOW RATE ON "),write(Link_in), nl,
        %Reset Flags when flow changes!
        Flag1 = 0, Flag2 = 0, Flag3 = 0, Flag4 = 0,
        nb_setval(albertSouthFlag,Flag1),
        nb_setval(callansSouthFlag,Flag2),
        nb_setval(heptonSouthFlag,Flag3),
        nb_setval(tradesNorthFlag,Flag4),

        do_plan_waiting(Rest,T, Rest2,Trig),
        !.
do_plan_waiting([R|Rest],T, [R|Rest2],Trig) :-
        do_plan_waiting(Rest,T, Rest2,Trig),
        !.
do_plan_waiting([],_,[],[]) :- !.


queue_monitor(InLink,FlagName,DYN,PLAN2, PLAN3) :-

     %Get the flag value
     nb_getval(FlagName,Flag),
     %Get the occupancy for the inlink
     member( equals( occupancy(InLink), Queue),  DYN),
     %If the flag is equal to 0, the stage lengths have not yet been changed
     ( Flag = 0 -> (

          %If the queue is more than 10, stage lengths need to be change
          ( Queue > 10 ->
              Delta is 3,
              %Call method to calculate and change stage lengths
              change_stage_lengths(InLink,_,_,Delta,PLAN3),
              %Stage lenghs have now been changed, so flag is updated so multiple               changes aren't made
              NewFlag is 1,
              nb_setval(FlagName,NewFlag)
              ;
              PLAN3 = PLAN2

              )
      );
        %If the flag isn't 0, the stage lengths have already changed so the plan can stay as it is
               PLAN3 = PLAN2
     ),
!.







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	simulation of events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

replace( L , X , Y , Z , R ) :-
  append(RowPfx,[Row|RowSfx],L),     % decompose the list-of-lists into a prefix, a list and a suffix
  length(RowPfx,X) ,                 % check the prefix length: do we have the desired list?
  append(ColPfx,[_|ColSfx],Row) ,    % decompose that row into a prefix, a column and a suffix
  length(ColPfx,Y) ,                 % check the prefix length: do we have the desired column?
  append(ColPfx,[Z|ColSfx],RowNew) , % if so, replace the column with its new value
  append(RowPfx,[RowNew|RowSfx],R)   % and assemble the transformed list-of-lists
  .


do_events :-
        assert(event_flag(0)),
        do_events1,
        retract(event_flag(N)),
        N > 0,
% keep repeating until no more events have happened
        do_events.
do_events.

do_events1 :-
        get_current_state_facts(N,T,ALL),
        event(_Name, _Types, Precons,Pos,Neg,Assign),
        satisfied1(Precons,ALL),
        satisfied2(Precons,ALL),
        write_out_if_its_a_trigger(T,Precons),
        retract(node(N,T,DYN1,PL,PR)),
        event_change(Pos,Neg,Assign,DYN1, DYN2),
        assert(node(N,T,DYN2,PL,PR)),
        retract(event_flag(Count)), N1 is Count+1, assert(event_flag(N1)),
        fail.
do_events1.

timed_var_clear(T):- (T > 90 ,T < 105) -> write('more than 100!'),!.

write_out_if_its_a_trigger(T,Precons) :-
        member(trigger(J),Precons),
        tell(user),nl,write("Time is: "),write(T), write(" seconds"),nl,
        write(J),write("  going to intergreen"),!.
write_out_if_its_a_trigger(_,_) :-!.

event_change(Pos,Neg,Assign,DYN1, DYN2) :-
        event_changeP(Pos, DYN1,DYN11),
        event_changeN(Neg, DYN11,DYN12),
        event_changeA(Assign, DYN12,DYN2),!.

event_changeP([],DYN,DYN) :- !.
event_changeP([A|B],DYN,[A | DYN2]) :-
        event_changeP(B,DYN,DYN2),!.

event_changeN([],DYN,DYN) :- !.
event_changeN([A|B],DYN,DYN2) :-
        delete(DYN,A, DYN1),
        event_changeN(B,DYN1,DYN2),!.

event_changeA([],DYN,DYN) :- !.
%Welcome to SWI-Prolog (threaded, 64 bits, version 8.0.3)
%SWI-Prolog comes with ABSOLUTELY NO WARRANTY. This is free software.
%Please run ?- license. for legal details.

%For online help and background, visit http://www.swi-prolog.0org
%For built-in help, use ?- help(Topic). or ?- apropos(Word).

event_changeA([assign(FUN,V)|B],DYN,[equals(FUN,V1)|DYN2]) :-
        get_statics(ST),
        append(DYN,ST,ALL),
        value_of(V,ALL,V1),
        delete(DYN,equals(FUN,_), DYN1),
        event_changeA(B,DYN1,DYN2),!.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	simulation of processes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_processes(Delta) :-
        get_current_state_facts(N1,T1,ALL),
        process(_Name,_Types,Pre,Ch),
% are the pres satisfied in the current state?
        satisfied1(Pre, ALL),
        satisfied2(Pre, ALL),
        retract(node(N1,T1,DYN1,PL,PR)),
% note that ALL keeps ALL the old state dynamic values whereas values in DYN1 .. DYN2.. are updated
        progress_change(Delta,DYN1,ALL,Ch, DYN2),
        assert(node(N1,T1,DYN2,PL,PR)),
       fail.
do_processes(_) :- !.

progress_change(Delta,DYN,ALL,[increase(FUN,EXP) | CR],  DYN2) :-
        value_of_op(Delta,ALL,EXP,V),
        increase_val(FUN,V,DYN,DYN1),
        progress_change(Delta,DYN1,ALL,CR,  DYN2),!.
progress_change(Delta,DYN,ALL,[decrease(FUN,EXP) | CR],  DYN2) :-
        value_of_op(Delta,ALL,EXP,V),
        decrease_val(FUN,V,DYN,DYN1),
	progress_change(Delta,DYN1,ALL,CR,  DYN2),!.
progress_change(_,DYN,_,[],  DYN) :- !.
progress_change(A,B,C,D,  E) :- nl,nl,nl,write('*** FAIL *** '),nl,
       write(A),nl,
       write(B),nl,
       write(C),nl,
       write(D),nl,
       write(E),nl,
       !.


increase_val(FUN,V,DYN,[equals(FUN,V2)|DYN1]) :-
        member(equals(FUN,V1),DYN),
        V2 is V1+V,
	delete(DYN,equals(FUN,_), DYN1),!.
decrease_val(FUN,V,DYN,[equals(FUN,V2)|DYN1]) :-
        member(equals(FUN,V1),DYN),
        V2 is V1-V,
	delete(DYN,equals(FUN,_), DYN1),!.


% deals with predicates first - leaves 'compare' till after
satisfied1([],_).
satisfied1( [compare(_, _, _) | R ], D) :-
        satisfied1(R,D).
satisfied1( [ X | R ], D) :-
        member(X,D),
        satisfied1(R, D).

% now deals with 'compare' - compare two expressions that will evaluate to a number
% get 2 and 3rd arg matched with equals(_, V and V1), then check V op V1.
satisfied2([],_).
satisfied2( [compare(Op, E1, E2) | R ], D) :-
        do_compare(Op,E1,E2, D),
        satisfied2(R,D).
satisfied2( [J | R ] , D) :-
        J \= compare(_,_,_),
        satisfied2(R, D).


do_compare(>,E1,E2, D) :-
         value_of(E1,D,V1),
         value_of(E2,D,V2),
         V1 > V2.
do_compare(<,E1,E2, D) :-
        value_of(E1,D,V1),
         value_of(E2,D,V2),
         V1 < V2.
do_compare(>=,E1,E2, D) :-
        value_of(E1,D,V1),
         value_of(E2,D,V2),
         V1 >= V2.
do_compare(=<,E1,E2, D) :-
        value_of(E1,D,V1),
         value_of(E2,D,V2),
         V1 =< V2.

% find the value of a function or constant - note if a function, some args could be vars
% this may need to be backtrackable - if we are evaluating in the preconditions it
% may depend on something like member(equals(greentime(X),ALL) where X needs to range through
value_of(op(-,E1,E2),ALL,E3) :-
      value_of(E1,ALL,V1),
      value_of(E2,ALL,V2),
      E3 is V1-V2.
value_of(op(+,E1,E2),ALL,E3) :-
      value_of(E1,ALL,V1),
      value_of(E2,ALL,V2),
      E3 is V1+V2.
value_of(op(*,E1,E2),ALL,E3) :-
      value_of(E1,ALL,V1),
      value_of(E2,ALL,V2),
      E3 is V1*V2.
value_of(op(/,E1,E2),ALL,E3) :-
      value_of(E1,ALL,V1),
      value_of(E2,ALL,V2),
      E3 is V1/V2.
value_of(E,_,E) :- float(E),!.
value_of(E,_,E) :- integer(E),!.
value_of(E,ALL,V) :- member(equals(E,V),ALL).

% deal with special case op(_,_,_)
% could put in value_of_op … X, X is integer or float, = X
value_of_op(Delta, _, op(*,t,X), V) :-
        integer(X),
        V is X*Delta,!.
value_of_op(Delta, _, op(*,t,X), V) :-
        float(X),
        V is X*Delta,!.
value_of_op(Delta,ALL,op(*,t,X), V) :-
        member(equals(X,VX),ALL),
        V is VX * Delta,!.
value_of_op(Delta,ALL,op(*,t,X), V) :-
        value_of(X,ALL,VX),
        V is VX * Delta,!.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% utilities

get_links(C) :- init(_,_,LC,_,_,_, _, _),member(link(C),LC),!.
get_stages(C) :- init(_,_,LC,_,_,_, _, _),member(stage(C),LC),!.
get_junctions(C) :- init(_,_,LC,_,_,_, _, _),member(junction(C),LC),!.
get_initial(II) :- init(_,_,_,J,I,_, _, _),append(J,I,II),!.
get_statics(STATICS) :-  statics(STATICS),!.
get_current_state_facts(N1,T1,ALL) :-
        get_statics(ST),
        node(N1,T1,DYN,_,_),
        append(ST,DYN,ALL),!.

memberNOC(X,[X|_]).
memberNOC(X,[_|R]) :- memberNOC(X,R).

%%%%% *************CHECKS   **************************************
/* check link one way */

check_links :- get_links(C), check(C),!.
check([ ]):- nl, write('All links checked for occupancy value'), nl,nl.
check([C1 | C]):- get_initial(I), member(equals(occupancy(C1),_), I), check(C), !.
check([C1 | C]):- write('error: link '),write(C1),write(' has no occupancy value'), nl, check(C).

check_links_max :- get_links(C), checkm(C),!.
checkm([ ]):- nl, write('All links checked for maximum occupancy value'), nl.
checkm([C1 | C]):- get_initial(I), member(equals(capacity(C1),_), I),checkm(C), !.
checkm([C1 | C]):- write('error: link '),write(C1),write(' has no maximum occupancy value'), nl, checkm(C).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% recording and writing out stuff .........

% write out each time step if trace is on:
log_trace(T,N,DYN,PL,PR) :-
         sim_trace(on),
         tracefile(Trace),
         tell(Trace),
         nl,nl,write("===================="),nl,
         write(N), write("   "), write(T),nl,
         get_statics(STAT),
% write trace records saturated links also......
         write_trace(N,T,DYN,STAT),
         write_traceX(DYN,STAT),
         write_traceY(DYN),
         nl,write("CURRENT PLAN:  "),nl,write(PL),
         nl,write("CURRENT PRESSES: "),nl,write(PR),
         tell(user),
         !.
log_trace(_,_,_,_,_,_) :-
         !.

write_end_stuff(STAT,PLAN,DYN,Sim_Time) :-
        nl,write(' END OF SIMULATION TIME' ),nl,
        member( equals( occupancy(station_east), Occ1),  DYN),
        member( equals( occupancy(station_west), Occ2),  DYN),
        member( equals( occupancy(stubbings_east), Occ3),  DYN),
        member( equals( occupancy(stubbings_west), Occ4),  DYN),
        member( equals( occupancy(albert_south), Occ5),  DYN),
        member( equals( occupancy(callans_south), Occ6),  DYN),
        member( equals( occupancy(trades_north), Occ7),  DYN),
        member( equals( occupancy(hepton_south), Occ8),  DYN),
        member(equals(turnrate(fake,outside,station_west),INrateE),STAT),
        member(equals(turnrate(fake,outside,stubbings_east),INrateW),STAT),
        member(equals(defaultgreentime(centre_stage0) , LongCentre),PLAN),
        member(equals(defaultgreentime(centre_stage1) , ShortCentre),PLAN),
        member(equals(defaultgreentime(vocation_stage0) , LongVocation),PLAN),
        member(equals(defaultgreentime(vocation_stage1) , ShortVocation),PLAN),
        member(equals(defaultgreentime(fox_stage0) , LongFox),PLAN),
        member(equals(defaultgreentime(fox_stage1) , ShortFox),PLAN),
        member(equals(defaultgreentime(pelican_stage0) , LongPelican),PLAN),
        member(equals(interlimit(pelican_stage0) , ShortPelican),STAT),
        domainmodel(DM),
        nl,write("  Domain Model is: "),write(DM),nl,
        initial(ISS),
        nl,write("  Initial State is: "),write(ISS),nl,
        nl,write("  Centre Timings: "),write(LongCentre),write(" / "),write(ShortCentre),
        write("  Vocation Timings: "),write(LongVocation),write(" / "),write(ShortVocation),
        write("  Fox Timings: "),write(LongFox),write(" / "),write(ShortFox),
        write("  Pelican Timings: "),write(LongPelican),write(" / "),write(ShortPelican),
        nl,nl,
        write("Saturation levels: "),nl,
        write_saturation_levels,nl,
        write("No of vehicles enter the region from the EAST: "),
        R is Sim_Time*INrateE - Occ2,write(R),nl,
        write("No of vehicles left the region in the WEST END: "),
        write(Occ4),nl,
        write("No of vehicles enter the region from the WEST: "),
        R1 is Sim_Time*INrateW - Occ3,write(R1),nl,
        write("No of vehicles left the region in the EAST END: "),
        write(Occ1),nl, nl,
        write("Queue waiting in the EAST: "),
        Occ2R is round(Occ2),write(Occ2R),nl,
        write("Queue waiting in the WEST: "),
        Occ3R is round(Occ3),write(Occ3R),nl,
        write("Queue waiting in ALBERT STREET: "),
        Occ5R is round(Occ5),write(Occ5R),nl,
        write("Queue waiting in HEPTON ROAD: "),
        Occ8R is round(Occ8),write(Occ8R),nl,
        write("Queue waiting in the CALLANS: "),
        Occ6R is round(Occ6),write(Occ6R),nl,
        write("Queue waiting in the TRADES: "),
        Occ7R is round(Occ7),write(Occ7R),nl,
        nl,
        !.

write_out_turns(STAT) :-
         member( equals(turnrate(A,B,C),TR), STAT),
         write(TR),write("  TURNRATE FOR: "),write(A),write("  "),write(B),write("  "),write(C),nl,
         fail.
write_out_turns(_) :- !.


assert_saturation_levels :-
        assert(saturation(park_east,0)),
        assert(saturation(park_west,0)),
        assert(saturation(new_east,0)),
        assert(saturation(new_west,0)),
        assert(saturation(market_east,0)),
        assert(saturation(market_west,0)).

write_saturation_levels :-
        saturation(X,LI),
        write("Time Saturated "), write(X), write(" is "),write(LI),nl,
        fail.
write_saturation_levels :- !.

write_list([ ]) :- nl,nl,write('END'),nl.
write_list([X|Y]) :- nl, write(X), write_list(Y),!.

write_listX([]) :- nl,write("["),nl,write("]"),nl.
write_listX([X]) :- nl,write("     "),write(X),nl,write("]"),nl.
write_listX([X|Y]) :- nl, write("     "),write(X), write(","),write_listX(Y),!.

write_trace(_,_,DYN,STAT) :-
         member(equals(occupancy(L),V), DYN),
         write_trace_1(L,V,STAT),
         fail.
write_trace(_,_,_,_).

write_traceX(DYN,ST) :-
         member( equals(greentime(J),GT), DYN),
         member(contains(J,S),ST), member(active(S),DYN),
         nl,write("GREEN TIME: "),write(J),write(" "),write("stage "),write(S),write(" "),write(GT),nl,
         fail.
write_traceX(_,_) :- !.
write_traceY(DYN) :-
         member( equals(queue(S),GT), DYN),
         nl,write("QUEUE: "),write(S),write(":   "),write(GT),nl,
         fail.
write_traceY(_) :- !.

record_sat(_,S) :- S < 85,!.
record_sat(L,_) :-
        retract(saturation(L,X)),
        X1 is X + 1,
        assert(saturation(L,X1)),!.

write_trace_1(L,V,STAT) :-
         member(equals(capacity(L),C), STAT),
         C < 9000,
         SAT is round(100*V/C),
%         V1 is round(V),
         V1 is V,
         record_sat(L,SAT),
         write(L),write("   OCCUPANCY:  "),write(V1), write(" SATURATION LEVEL:   "),write(SAT),nl,!.
write_trace_1(L,V,STAT) :-
         member(equals(capacity(L),C), STAT),
         C >= 9000,
%         V1 is round(V),
         V1 is V,
         write(L),write("   OCCUPANCY:  "),write(V1),nl,!.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
