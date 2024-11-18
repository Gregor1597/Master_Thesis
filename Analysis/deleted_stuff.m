%%
diff_C_FM1 = CoM(1,:)- foot_middle(1,:);
diff_C_FM2 = CoM(2,:)- foot_middle(2,:);
euclid = sqrt(diff_C_FM1.^2 + diff_C_FM2.^2);
first = true;
trials = [];
if (first)
    walking = euclid > 15;
    start = find(walking);
    disp(start(1,1));
   % e = CoM(1,start(1,1):end);
    %ends = CoM(1,:) < CoM(1, start(1,1));
    ends = abs(diff(CoM(1,:))) < 0.00095;
    e = find(ends);
    disp(e(start(1,1)+1));
    ttt = CoM(1,start(1,1):e(start(1,1)+1));
    trials{1} = ttt;
    first = false;
%else?
    cut_e = euclid(1, e(start(1,1)+1):end) > 15;
    cut_c = CoM(:, e(start(1,1)+1):end);
    start2 = find(cut_e);
    disp(start2(1,1));
    %ends2 = cut_c(1,:) < cut_c(1, start2(1,1));
    ends2 = abs(diff(cut_c(1,:))) < 0.00095;
    e2 = find(ends2);
    disp(e2(start2(1,1)+1));
    ttt2 = cut_c(1, start2(1,1):e2(start2(1,1)+1));
    trials{2} = ttt2;
    
    cut_e3 = euclid(1, e2(start2(1,1)+1):end) > 15;
    cut_c3 = CoM(:, e2(start2(1,1)+1):end);
    start3 = find(cut_e3);
    disp(start3(1,1));
    %ends2 = cut_c(1,:) < cut_c(1, start2(1,1));
    ends3 = abs(diff(cut_c3(1,:))) < 0.00095;
    e3 = find(ends3);
    disp(e3(start3(1,1)+1));
    ttt3 = cut_c3(1, start3(1,1):e3(start3(1,1)+1));
    trials{3} = ttt3;
end

%%
%% test trails schneiden
ms = CoM(1,:);
Df = diff(ms);
plot(ms)
hold on
plot(Df*200)
%%
plot(CoM(1,:),CoM(2,:))
hold on 
plot(foot_middle(1,:), foot_middle(2,:))