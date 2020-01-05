
% initial state

init(
% problem name:
hebden_bridge,

% linked with domain model:
urbantraffic,

% list of type_name( list of objects of type_name ) 
[
     junction([fox,vocation,centre,pelican]),
     link([new_west,market_east,market_west,callans_south,callans_north,park_west,new_east,station_west,station_east,stubbings_west,stubbings_east,albert_north,albert_south,trades_south, hepton_north, hepton_south, trades_north,park_east,outside]),
     stage([pelican_stage0,fox_stage0,fox_stage1,centre_stage0,centre_stage1,vocation_stage0,vocation_stage1, fake]),
     section([new_s1,new_s2,new_s3,new_s4,new_s5,market_s1,market_s2,market_s3,market_s4,market_s5,market_s6,market_s7,market_s8,market_s9,market_s10,park_s1,park_s2,park_s3,park_s4,park_s5,new_east_s1,new_east_s2,new_east_s3,new_east_s4,new_east_s5,market_east_s1,market_east_s2,market_east_s3,market_east_s4,market_east_s5,market_east_s6,market_east_s7,market_east_s8,market_east_s9,market_east_s10,park_east_s1,park_east_s2,park_east_s3,park_east_s4,park_east_s5])
],

%  dynamic facts - ones that can be changed by the process simulation 
[
% is true if he stage is on green
% 'fake' is a 'stage' that is always active - like traffic flowing into the boundary
     active(fake),
     active(pelican_stage0),
     active(centre_stage0),
     active(fox_stage0),
     active(vocation_stage0),

% this holds the number of PCUs at the end of the link 
     equals(queue(new_west),0.0),
     equals(queue(park_west),0.0),
     equals(queue(market_west),0.0),
     equals(queue(new_east),0.0),
     equals(queue(park_east),0.0),
     equals(queue(market_east),0.0),

% sections are used to regulate flow through a link - this is a pointer to a section
     current(new_s1,new_west),
     current(park_s1,park_west),
     current(market_s1,market_west),
     current(new_east_s1,new_east),
     current(park_east_s1,park_east),
     current(market_east_s1,market_east),

% this records the amount of time traffic in a section have been waiting in link 
     equals(sectiontime(new_west),0),
     equals(sectiontime(park_west),0),
     equals(sectiontime(market_west),0),
     equals(sectiontime(new_east),0),
     equals(sectiontime(park_east),0),
     equals(sectiontime(market_east),0),

% amount of time a light has been on green
     equals(greentime(fox),0),
     equals(intertime(fox),0),
     equals(greentime(centre),0),
     equals(intertime(centre),0),
     equals(greentime(pelican),0),
     equals(intertime(pelican),0),
     equals(greentime(vocation),0),
     equals(intertime(vocation),0),

% amount of PCUs in a link
     equals(occupancy(park_east),0.0),
     equals(occupancy(park_west),0.0),
     equals(occupancy(new_east),0.0),
     equals(occupancy(new_west),0.0),
     equals(occupancy(market_east),0.0),
     equals(occupancy(market_west),0.0),
     equals(occupancy(stubbings_east),0.0),
     equals(occupancy(stubbings_west),0.0),
     equals(occupancy(station_east),0.0),
     equals(occupancy(station_west),0.0),
     equals(occupancy(callans_south),0.0),
     equals(occupancy(callans_north),0.0),
     equals(occupancy(hepton_south),0.0),
     equals(occupancy(hepton_north),0.0),
     equals(occupancy(albert_north),0.0),
     equals(occupancy(albert_south),0.0),
     equals(occupancy(trades_north),0.0),
     equals(occupancy(trades_south),0.0),
     equals(occupancy(outside),99999.0),

% amount of PCUs in a section
     equals(store(new_s1),0.0),
     equals(store(new_s2),0.0),
     equals(store(new_s3),0.0),
     equals(store(new_s4),0.0),
     equals(store(new_s5),0.0),
     equals(store(market_s1),0.0),
     equals(store(market_s2),0.0),
     equals(store(market_s3),0.0),
     equals(store(market_s4),0.0),
     equals(store(market_s5),0.0),
     equals(store(market_s6),0.0),
     equals(store(market_s7),0.0),
     equals(store(market_s8),0.0),
     equals(store(market_s9),0.0),
     equals(store(market_s10),0.0),
     equals(store(park_s1),0.0),
     equals(store(park_s2),0.0),
     equals(store(park_s3),0.0),
     equals(store(park_s4),0.0),
     equals(store(park_s5),0.0),
     equals(store(new_east_s1),0.0),
     equals(store(new_east_s2),0.0),
     equals(store(new_east_s3),0.0),
     equals(store(new_east_s4),0.0),
     equals(store(new_east_s5),0.0),
     equals(store(market_east_s1),0.0),
     equals(store(market_east_s2),0.0),
     equals(store(market_east_s3),0.0),
     equals(store(market_east_s4),0.0),
     equals(store(market_east_s5),0.0),
     equals(store(market_east_s6),0.0),
     equals(store(market_east_s7),0.0),
     equals(store(market_east_s8),0.0),
     equals(store(market_east_s9),0.0),
     equals(store(market_east_s10),0.0),
     equals(store(park_east_s1),0.0),
     equals(store(park_east_s2),0.0),
     equals(store(park_east_s3),0.0),
     equals(store(park_east_s4),0.0),
     equals(store(park_east_s5),0.0)
],

% list of STATIC facts - ones that are not changed by the simulation process
[

% length of intergreen after a stage
     equals(interlimit(centre_stage0),5),
     equals(interlimit(centre_stage1),5),
     equals(interlimit(vocation_stage0),5),
     equals(interlimit(vocation_stage1),5),
     equals(interlimit(fox_stage0),5),
     equals(interlimit(fox_stage1),5),
     equals(interlimit(pelican_stage0),15),

% connects junctions with stages
     contains(centre,centre_stage0),
     contains(centre,centre_stage1),
     contains(fox,fox_stage0),
     contains(fox,fox_stage1),
     contains(vocation,vocation_stage0),
     contains(vocation,vocation_stage1),
     contains(pelican,pelican_stage0),

% maximum capacity of a link in PCUs
     equals(capacity(park_west),20.0),
     equals(capacity(park_east),20.0),
     equals(capacity(new_east),20.0),
     equals(capacity(new_west),20.0),
     equals(capacity(market_east),40.0),
     equals(capacity(market_west),40.0),
     equals(capacity(callans_north),9999.0),
     equals(capacity(callans_south),9999.0),
     equals(capacity(stubbings_west),9999.0),
     equals(capacity(stubbings_east),9999.0),
     equals(capacity(station_east),9999.0),
     equals(capacity(station_west),9999.0),
     equals(capacity(albert_north),9999.0),
     equals(capacity(albert_south),9999.0),
     equals(capacity(trades_north),9999.0),
     equals(capacity(trades_south),9999.0),
     equals(capacity(hepton_north),9999.0),
     equals(capacity(hepton_south),9999.0),
     equals(capacity(outside),99999.0),

% this holds the value of the flow rate during a stage of a junction 
% between two links 
% this first specifies flow from outside

     equals(turnrate(fake,outside,hepton_south),0.07),
     equals(turnrate(fake,outside,albert_south),0.05),
     equals(turnrate(fake,outside,callans_south),0.05),
     equals(turnrate(fake,outside,trades_north),0.05),
     equals(turnrate(fake,outside,stubbings_east),0.25),
     equals(turnrate(fake,outside,station_west),0.25),

     equals(turnrate(fox_stage0,stubbings_east,market_east),0.4),
     equals(turnrate(fox_stage0,market_west,stubbings_west),0.4),
     equals(turnrate(fox_stage0,stubbings_east,hepton_north),0.01),
     equals(turnrate(fox_stage0,market_west,hepton_north),0.01),
     equals(turnrate(fox_stage1,hepton_south, market_east),0.1),
     equals(turnrate(fox_stage1,hepton_south, stubbings_west),0.1),



     equals(turnrate(vocation_stage0,station_west,new_west),0.4),
     equals(turnrate(vocation_stage0,new_east,station_east),0.4),
     equals(turnrate(vocation_stage0,station_west,albert_north),0.03),
     equals(turnrate(vocation_stage0,new_east,albert_north),0.03),
     equals(turnrate(vocation_stage1,albert_south,new_west),0.2),
     equals(turnrate(vocation_stage1,albert_south,station_east),0.2),

     equals(turnrate(centre_stage0,market_east,callans_north),0.015),
     equals(turnrate(centre_stage0,market_east,trades_south),0.015),
     equals(turnrate(centre_stage0,market_east,park_east),0.4),
     equals(turnrate(centre_stage0,park_west,market_west),0.4),
     equals(turnrate(centre_stage0,park_west,callans_north),0.015),
     equals(turnrate(centre_stage0,park_west,trades_south),0.015),

     equals(turnrate(centre_stage1,callans_south,park_east),0.15),
     equals(turnrate(centre_stage1,callans_south,market_west),0.05),
     equals(turnrate(centre_stage1,callans_south,trades_south),0.15),
     equals(turnrate(centre_stage1,trades_north,park_east),0.05),
     equals(turnrate(centre_stage1,trades_north,market_west),0.15),
     equals(turnrate(centre_stage1,trades_north,callans_north),0.15),

     equals(turnrate(pelican_stage0,park_east,new_east),0.4),
     equals(turnrate(pelican_stage0,new_west,park_west),0.4),
%
     equals(mingreentime(centre_stage0),10),
     equals(maxgreentime(centre_stage0),60),
     equals(mingreentime(centre_stage1),5),
     equals(maxgreentime(centre_stage1),20),
     equals(mingreentime(pelican_stage0),20),
     equals(maxgreentime(pelican_stage0),100000),
     equals(mingreentime(fox_stage0),20),
     equals(maxgreentime(fox_stage0),60),
     equals(mingreentime(fox_stage1),5),
     equals(maxgreentime(fox_stage1),20),
     equals(mingreentime(vocation_stage0),10),
     equals(maxgreentime(vocation_stage0),60),
     equals(mingreentime(vocation_stage1),5),
     equals(maxgreentime(vocation_stage1),25),
% specifies order of stage
     next(pelican_stage0,pelican_stage0),
     next(centre_stage0,centre_stage1),
     next(centre_stage1,centre_stage0),
     next(vocation_stage0,vocation_stage1),
     next(vocation_stage1,vocation_stage0),
     next(fox_stage0,fox_stage1),
     next(fox_stage1,fox_stage0),
% specifies order of sections
     nextsection(new_s1,new_s2),
     nextsection(new_s2,new_s3),
     nextsection(new_s3,new_s4),
     nextsection(new_s4,new_s5),
     nextsection(new_s5,new_s1),
     nextsection(market_s1,market_s2),
     nextsection(market_s2,market_s3),
     nextsection(market_s3,market_s4),
     nextsection(market_s4,market_s5),
     nextsection(market_s5,market_s6),
     nextsection(market_s6,market_s7),
     nextsection(market_s7,market_s8),
     nextsection(market_s8,market_s9),
     nextsection(market_s9,market_s10),
     nextsection(market_s10,market_s1),
     nextsection(park_s1,park_s2),
     nextsection(park_s2,park_s3),
     nextsection(park_s3,park_s4),
     nextsection(park_s4,park_s5),
     nextsection(park_s5,park_s1),
     nextsection(new_east_s1,new_east_s2),
     nextsection(new_east_s2,new_east_s3),
     nextsection(new_east_s3,new_east_s4),
     nextsection(new_east_s4,new_east_s5),
     nextsection(new_east_s5,new_east_s1),
     nextsection(market_east_s1,market_east_s2),
     nextsection(market_east_s2,market_east_s3),
     nextsection(market_east_s3,market_east_s4),
     nextsection(market_east_s4,market_east_s5),
     nextsection(market_east_s5,market_east_s6),
     nextsection(market_east_s6,market_east_s7),
     nextsection(market_east_s7,market_east_s8),
     nextsection(market_east_s8,market_east_s9),
     nextsection(market_east_s9,market_east_s10),
     nextsection(market_east_s10,market_east_s1),
     nextsection(park_east_s1,park_east_s2),
     nextsection(park_east_s2,park_east_s3),
     nextsection(park_east_s3,park_east_s4),
     nextsection(park_east_s4,park_east_s5),
     nextsection(park_east_s5,park_east_s1),
% this records the length of time traffic stay in a section  in link I 
    equals(maxsectiontime(new_west),3),
     equals(maxsectiontime(park_west),3),
     equals(maxsectiontime(market_west),5),
     equals(maxsectiontime(new_east),3),
     equals(maxsectiontime(park_east),3),
     equals(maxsectiontime(market_east),5)
   
],

[
% these form the Default Plan !
% they give the green time for each stage in the 4 traffic signals
     equals(defaultgreentime(centre_stage0),30),
     equals(defaultgreentime(centre_stage1),20),
     equals(defaultgreentime(vocation_stage0),20),
     equals(defaultgreentime(vocation_stage1),20),
     equals(defaultgreentime(fox_stage0),20),
     equals(defaultgreentime(fox_stage1),20),
     equals(defaultgreentime(pelican_stage0),100000),
     change_flow(200,trades_north,0.12),
     change_flow(600,hepton_south,0.18),
     change_flow(400,callans_south,0.15),

     wait(pelican, 25)
],

% button presses with junction (pelican) and time pressed:
[
press(2,pelican), press(3,pelican), press(14,pelican), press(100,pelican),
press(150,pelican), press(200,pelican), press(270,pelican), press(310,pelican),
press(335,pelican), press(365,pelican), press(390,pelican), press(450,pelican),
press(500,pelican), press(600,pelican), press(650,pelican), press(700,pelican),
press(750,pelican), press(800,pelican), press(830,pelican), press(860,pelican),
press(900,pelican), press(1000,pelican), press(1100,pelican)
],

% goal - list of predicates - only relevant if used with a plan generation program
[
     compare(<,occupancy(station_west),200.0)
]

).

