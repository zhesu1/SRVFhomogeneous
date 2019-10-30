function [v, qv] = S2_get_initial_velocity(f1,q1,f2,q2)

v = logm(f1\f2);
qv = q2-q1;