% Tests on Artificial Networks
%SBM graphs
function test_ART

%informative case
for i=[2,2.3,2.5,2.8,3] %p_in/p_out
    run(0.1/i,2) %k=2 layers
    run(0.1/i,3) %k=3 layers
end

%noisy case
for i=[2.5,2.8,3,3.3,3.5]
    run_n(0.1/i,2,1) %2 info layers + 1 noisy layer
    run_n(0.1/i,2,2) %2 info layers + 2 noisy layers 
end
%}
end
