% Tests on Artificial Networks
%LFR graphs
function test_ART_LFR

%informative case
run_LFR(2) %k=2 layers
run_LFR(3) %k=3 layers

%noisy case
run_n_LFR(2,1) %2 info layers + 1 noisy layer
run_n_LFR(2,2) %2 info layers + 2 noisy layers 
%}
end