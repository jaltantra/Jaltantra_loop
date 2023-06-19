Sets
	nodes /1,	2,	3,	4,	5,	6,	7 /
	pipes /0,	1,	2,	3,	4,	5,	6,	7,	8,	9,	10,	11,	12,	13 /
	src(nodes) /1/;
alias (src,srcs);
alias (nodes,j) ;
Set arcs(nodes,j) //;
Set F_arcs(nodes,j) /1.2,  2.4,  2.3,  3.5,  4.5,  4.6,  7.5,  6.7/;

Parameters
	Len(nodes,j) //
	F_L(nodes,j) /1	.2	1000.0,  2	.4	1000.0,  2	.3	1000.0,  3	.5	1000.0,  4	.5	1000.0,  4	.6	1000.0,  7	.5	1000.0,  6	.7	1000.0/
	E(nodes) /1  210.0, 2  150.0, 3  160.0, 4  155.0, 5  150.0, 6  165.0, 7  160.0/
	P(nodes) /1  0.0, 2  30.0, 3  30.0, 4  30.0, 5  30.0, 6  30.0, 7  30.0/
	D(nodes) /1  -311.1087, 2  27.7777, 3  27.777, 4  33.333, 5  75.0, 6  91.666, 7  55.555/
	dia(pipes) /0  25.4, 1  50.8, 2  76.2, 3  101.6, 4  152.4, 5  203.2, 6  254.0, 7  304.8, 8  355.6, 9  406.4, 10  457.2, 11  508.0, 12  558.8, 13  609.6/
	C(pipes) /0  2.0, 1  5.0, 2  8.0, 3  11.0, 4  16.0, 5  23.0, 6  32.0, 7  50.0, 8  60.0, 9  90.0, 10  130.0, 11  170.0, 12  300.0, 13  550.0/
	R(pipes) /0  130.0, 1  130.0, 2  130.0, 3  130.0, 4  130.0, 5  130.0, 6  130.0, 7  130.0, 8  130.0, 9  130.0, 10  130.0, 11  130.0, 12  130.0, 13  130.0/
	F_d(nodes,j) /1	.2	600.0,  2	.4	600.0,  2	.3	600.0,  3	.5	600.0,  4	.5	600.0,  4	.6	600.0,  7	.5	600.0,  6	.7	600.0/
	F_R(nodes,j) /1	.2	130.0,  2	.4	130.0,  2	.3	130.0,  3	.5	130.0,  4	.5	130.0,  4	.6	130.0,  7	.5	130.0,  6	.7	130.0/;

Scalar omega  /10.68/;
Scalar bnd ;
Scalar qm;
Scalar q_M;

bnd = sum(src,D(src));
q_M=-bnd;
qm=0;

Variable l(nodes,j,pipes); 
l.lo(nodes,j,pipes)= 0;

Variable q1(nodes,j);
q1.lo(nodes,j)=qm;
q1.up(nodes,j)=q_M;

Variable q2(nodes,j);
q2.lo(nodes,j)=qm;
q2.up(nodes,j)=q_M;

Variables z;

Variable h(nodes);

Equations cost "objective function",bound1(nodes,j,pipes),cons1(nodes),cons2(nodes),cons3(nodes,j),cons5(src), cons4(nodes,j), cons6(nodes,j) ;
cost..  z=e=sum(arcs(nodes,j),sum(pipes,l(arcs,pipes)*c(pipes)));
bound1(nodes,j,pipes)$arcs(nodes,j).. l(nodes,j,pipes) =l= Len(nodes,j);
cons1(nodes).. sum(arcs(j,nodes),(q1(arcs)-q2(arcs))) =e= sum(arcs(nodes,j),(q1(arcs)-q2(arcs))) + D(nodes);
cons2(nodes).. h(nodes) =g= E(nodes) + P(nodes);
cons3(arcs(nodes,j)).. h(nodes)-h(j)=e=sum(pipes,(((q1(arcs)*0.001)**1.852 - (q2(arcs)*0.001)**1.852)*omega*l(arcs,pipes))/((R(pipes)**1.852)*(dia(pipes)/1000)**4.87));
cons4(arcs(nodes,j)).. sum(pipes,l(arcs,pipes)) =e=Len(arcs);
cons5(src)..  h(src)=e= sum(srcs,E(srcs));
cons6(arcs(nodes,j)).. q1(arcs)*q2(arcs) =l= q_M*qm;

model m2  /all/  ;
solve m2 using minlp minimizing z ;
