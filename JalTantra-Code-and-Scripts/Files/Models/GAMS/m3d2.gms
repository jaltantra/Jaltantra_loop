Sets
    nodes   /Node1, Node2, Node3, Node4, Node5, Node6, Node7, Node8, Node9, Node10, Node11, Node12, Node13, Node14, Node15, Node16,
             Node17, Node18, Node19, Node20, Node21, Node22, Node23, Node24, Node25, Node26, Node27, Node28, Node29, Node30, Node31, Node32/
    cycle   /cycle1, cycle2, cycle3/
    links   /Pipe_1_2, Pipe_2_3, Pipe_3_4, Pipe_4_5, Pipe_5_6, Pipe_6_7, Pipe_7_8, Pipe_8_9, Pipe_9_10, Pipe_10_11, Pipe_11_12, Pipe_12_13,
             Pipe_10_14, Pipe_14_15, Pipe_15_16, Pipe_17_16, Pipe_18_17, Pipe_19_18, Pipe_3_19, Pipe_3_20, Pipe_20_21, Pipe_21_22, Pipe_20_23,
             Pipe_23_24, Pipe_24_25, Pipe_26_25, Pipe_27_26, Pipe_16_27, Pipe_23_28, Pipe_28_29, Pipe_29_30, Pipe_30_31, Pipe_32_31, Pipe_25_32/
    pipes   /1, 2, 3, 4, 5, 6/
    src(nodes) /Node1/;

Parameters
    link_length(links) / Pipe_1_2   100, Pipe_2_3   1350, Pipe_3_4   900, Pipe_4_5   1150, Pipe_5_6   1450, Pipe_6_7   450, Pipe_7_8   850, Pipe_8_9   850, Pipe_9_10   800,
                        Pipe_10_11   950, Pipe_11_12   1200, Pipe_12_13   3500, Pipe_10_14   800, Pipe_14_15   500, Pipe_15_16   550, Pipe_17_16   2730, Pipe_18_17   1750,
                        Pipe_19_18   800, Pipe_3_19   400, Pipe_3_20   2200, Pipe_20_21   1500, Pipe_21_22   500, Pipe_20_23   2650, Pipe_23_24   1230, Pipe_24_25   1300,
                        Pipe_26_25   850, Pipe_27_26   300, Pipe_16_27   750, Pipe_23_28   1500, Pipe_28_29   2000, Pipe_29_30   1600, Pipe_30_31   150, Pipe_32_31   860, Pipe_25_32   950/
    elevation(nodes) /Node1   100, Node2   0, Node3   0, Node4   0, Node5   0, Node6   0, Node7   0, Node8   0, Node9   0, Node10   0, Node11   0, Node12   0,
                    Node13   0, Node14   0, Node15   0, Node16   0, Node17   0, Node18   0, Node19   0, Node20   0, Node21   0, Node22   0, Node23   0, 
                    Node24   0, Node25   0, Node26   0, Node27   0, Node28   0, Node29   0, Node30   0, Node31   0, Node32   0/
    pressure(nodes) /Node1   0, Node2   30, Node3   30, Node4   30, Node5   30, Node6   30, Node7   30, Node8   30, Node9   30, Node10   30, Node11   30, Node12   30,
                Node13   30, Node14   30, Node15   30, Node16   30, Node17   30, Node18   30, Node19   30, Node20   30, Node21   30, Node22   30, Node23   30, 
                Node24   30, Node25   30, Node26   30, Node27   30, Node28   30, Node29   30, Node30   30, Node31   30, Node32   30/
    demand(nodes) /Node1   -5538.8658000000005, Node2   247.222, Node3   236.111, Node4   36.111, Node5   201.388, Node6   279.1666, Node7   375, Node8   152.778,
                Node9   145.833, Node10   145.833, Node11   138.88, Node12   155.55, Node13   261.11, Node14   170.833, Node15   77.777, Node16   86.111,
                Node17   240.277, Node18   373.611, Node19   16.666, Node20   354.1666, Node21   258.333, Node22   134.7222, Node23   290.2777, 
                Node24   227.777, Node25   47.222, Node26   250, Node27   102.777, Node28   80.555, Node29   100, Node30   100, Node31   29.1666, Node32   223.6111/
    diameter(pipes) /1   304.8, 2   406.4, 3   508, 4   609.6, 5   762, 6   1016/
    Cost(pipes) /1   45.72617132, 2   70.4, 3   98.387, 4   129.33, 5   180.748, 6   278.28/
    Roughness(pipes) /1   130, 2   130, 3   130, 4   130, 5   130, 6   130/
    sourcehead(src) /Node1 100/;

Table F(nodes,links) "Flow Direction Matrix"
                Pipe_1_2 Pipe_2_3 Pipe_3_4 Pipe_4_5 Pipe_5_6 Pipe_6_7 Pipe_7_8 Pipe_8_9 Pipe_9_10 Pipe_10_11 Pipe_11_12 Pipe_12_13 Pipe_10_14 Pipe_14_15 Pipe_15_16 Pipe_17_16 Pipe_18_17 Pipe_19_18 Pipe_3_19 Pipe_3_20 Pipe_20_21 Pipe_21_22 Pipe_20_23 Pipe_23_24 Pipe_24_25 Pipe_26_25 Pipe_27_26 Pipe_16_27 Pipe_23_28 Pipe_28_29 Pipe_29_30 Pipe_30_31 Pipe_32_31 Pipe_25_32
        Node1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node2   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node3   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node4   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node5   0   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node6   0   0   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node7   0   0   0   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node8   0   0   0   0   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node9   0   0   0   0   0   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node10   0   0   0   0   0   0   0   0   1   -1   0   0   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node11   0   0   0   0   0   0   0   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node12   0   0   0   0   0   0   0   0   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node13   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node14   0   0   0   0   0   0   0   0   0   0   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node15   0   0   0   0   0   0   0   0   0   0   0   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node16   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   -1   0   0   0   0   0   0
        Node17   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node18   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node19   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node20   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   -1   0   -1   0   0   0   0   0   0   0   0   0   0   0
        Node21   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   -1   0   0   0   0   0   0   0   0   0   0   0   0
        Node22   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0
        Node23   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   -1   0   0   0   0   -1   0   0   0   0   0
        Node24   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   -1   0   0   0   0   0   0   0   0   0
        Node25   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   -1
        Node26   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1   1   0   0   0   0   0   0   0
        Node27   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1   1   0   0   0   0   0   0
        Node28   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   -1   0   0   0   0
        Node29   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   -1   0   0   0
        Node30   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   -1   0   0
        Node31   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0
        Node32   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1   1;


Table S(nodes,links) "Matrix for flow Direction in Spanning Tree"
                Pipe_1_2 Pipe_2_3 Pipe_3_4 Pipe_4_5 Pipe_5_6 Pipe_6_7 Pipe_7_8 Pipe_8_9 Pipe_9_10 Pipe_10_11 Pipe_11_12 Pipe_12_13 Pipe_10_14 Pipe_14_15 Pipe_15_16 Pipe_17_16 Pipe_18_17 Pipe_19_18 Pipe_3_19 Pipe_3_20 Pipe_20_21 Pipe_21_22 Pipe_20_23 Pipe_23_24 Pipe_24_25 Pipe_26_25 Pipe_27_26 Pipe_16_27 Pipe_23_28 Pipe_28_29 Pipe_29_30 Pipe_30_31 Pipe_32_31 Pipe_25_32
        Node1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node2   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node3   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node4   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node5   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node6   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node7   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node8   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node9   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node10   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node11   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node12   1   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node13   1   1   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node14   1   1   1   1   1   1   1   1   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node15   1   1   1   1   1   1   1   1   1   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node16   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node17   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node18   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node19   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node20   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node21   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0
        Node22   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0
        Node23   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   -1   -1   1   1   1   0   0   0   0   0   0
        Node24   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   -1   1   1   1   0   0   0   0   0   0
        Node25   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   0   0   0   0
        Node26   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0
        Node27   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0
        Node28   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   -1   -1   1   1   1   1   0   0   0   0   0
        Node29   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   -1   -1   1   1
        Node30   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   0   -1   1   1
        Node31   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   0   0   1   1
        Node32   1   1   1   1   1   1   1   1   1   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   0   0   0   1;


Table C(cycle,links) "Cycle Flow Direction Matrix"
                Pipe_1_2 Pipe_2_3 Pipe_3_4 Pipe_4_5 Pipe_5_6 Pipe_6_7 Pipe_7_8 Pipe_8_9 Pipe_9_10 Pipe_10_11 Pipe_11_12 Pipe_12_13 Pipe_10_14 Pipe_14_15 Pipe_15_16 Pipe_17_16 Pipe_18_17 Pipe_19_18 Pipe_3_19 Pipe_3_20 Pipe_20_21 Pipe_21_22 Pipe_20_23 Pipe_23_24 Pipe_24_25 Pipe_26_25 Pipe_27_26 Pipe_16_27 Pipe_23_28 Pipe_28_29 Pipe_29_30 Pipe_30_31 Pipe_32_31 Pipe_25_32
        cycle1   0   0   -1   -1   -1   -1   -1   -1   -1   0   0   0   -1   -1   -1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
        cycle2   0   0   -1   -1   -1   -1   -1   -1   -1   0   0   0   -1   -1   -1   0   0   0   0   1   0   0   1   1   1   -1   -1   -1   0   0   0   0   0   0
        cycle3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1   -1   0   0   0   1   1   1   1   -1   -1;


Scalar omega  /10.68/;
Scalar bnd ;
Scalar qm;
Scalar q_M;

bnd = sum(src,demand(src));
q_M=-bnd;
qm=10**(-1);

Variable l(links,pipes);
l.lo(links,pipes) = 0;

Variable q(links);
q.lo(links) = -q_M;
q.up(links) = q_M;

Variables z;

Equations obj "objective function",bound1(links,pipes),cons1(nodes),cons2(cycle),cons3a(nodes),cons3b(nodes), cons4(links), cons5(links);
obj.. z=e=sum(links,sum(pipes,l(links,pipes)*Cost(pipes)));

bound1(links,pipes).. l(links,pipes) =l= link_length(links);

cons1(nodes).. sum(links,F(nodes,links)*q(links)) =e= demand(nodes);

cons2(cycle).. sum(links,sum(pipes,C(cycle,links)*omega*l(links,pipes)*(q(links)*0.001)*(abs((q(links)*0.001))**0.852)/((Roughness(pipes)**1.852)*((diameter(pipes)/1000)**4.87)))) =e= 0;

cons3a(nodes).. sum(links,sum(pipes,S(nodes,links)*omega*l(links,pipes)*(q(links)*0.001)*(abs((q(links)*0.001))**0.852)/((Roughness(pipes)**1.852)*((diameter(pipes)/1000)**4.87)))) =l= sum(src,sourcehead(src))-elevation(nodes)-pressure(nodes);
cons3b(nodes).. sum(links,sum(pipes,S(nodes,links)*omega*l(links,pipes)*(q(links)*0.001)*(abs((q(links)*0.001))**0.852)/((Roughness(pipes)**1.852)*((diameter(pipes)/1000)**4.87)))) =g= 0;

cons4(links).. qm =l= abs(q(links));
cons5(links).. sum(pipes,l(links,pipes)) =e= link_length(links);

model m3d1  /all/  ;
solve m3d1 using dnlp minimizing z ;

