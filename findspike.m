function spike = findspike(V1) 
    spike = V1(2:end-1)>V1(1:end-2)&V1(2:end-1)>V1(3:end)&V1(2:end-1)>0.01;
    spike = [false spike false];
end

