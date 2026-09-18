
% 棒球軌跡【三變數聯合鑑別】系統 (包含隱藏仰角偵測)
clear; clc;

%% 1. 輸入全部實驗數據 (Ground Truth)
% 【數據A：高轉速，17m內彈一次】
v0_actual_High = 18.67;      % 高轉速組實測初速 (m/s)
dist_H_1st_actual = 12.48;   % 高轉速組實測第一落點 (公尺)
% 【數據B：低轉速，17m內彈兩次】
v0_actual_Low = 13.14;       % 低轉速組實測初速 (m/s)
dist_L_1st_actual = 8.61;    % 低轉速組實測第一落點 (公尺)
dist_L_2nd_actual = 13.77;   % 低轉速組實測第二落點 (公尺)
fprintf('=== 全實驗數據匯入 ===\n');
fprintf('[高轉速] 實測初速: %.2f m/s, 第一落點: %.2f m\n', v0_actual_High, dist_H_1st_actual);
fprintf('[低轉速] 實測初速: %.2f m/s, 第一落點: %.2f m, 第二落點: %.2f m\n', v0_actual_Low, dist_L_1st_actual, dist_L_2nd_actual);
%% 2. 第一階段：聯合鑑別空氣阻力 (Cd) 與 隱藏仰角 (Theta)
fprintf('\n=== [階段一] 聯合鑑別空氣阻力與隱藏仰角 ===\n');
% x0_phase1(1) = 空氣阻力 Cd (猜測 0.3)
% x0_phase1(2) = 發射仰角 Theta 角度 (猜測 5 度)
x0_phase1 = [0.3, 5.0]; 
options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);
[opt_params_1, state_at_bounce_Low] = fminsearch_Phase1(x0_phase1, v0_actual_High, dist_H_1st_actual, v0_actual_Low, dist_L_1st_actual, options);
optimized_Cd = opt_params_1(1);
optimized_theta = opt_params_1(2);
fprintf('-> 最佳空氣阻力係數 (Cd): %.4f\n', optimized_Cd);
fprintf('-> 偵測到機台隱藏仰角:   %.2f 度\n', optimized_theta);
fprintf('-> (低轉速組) 落地瞬間水平速度: %.2f m/s, 垂直速度: %.2f m/s\n', state_at_bounce_Low(3), state_at_bounce_Low(4));
%% 3. 第二階段：鑑別草地恢復係數 (e_y, e_x)
fprintf('\n=== [階段二] 鑑別草地恢復係數 ===\n');
% x0_e(1) = 垂直恢復係數 e_y (猜測 0.4)
% x0_e(2) = 水平摩擦衰減 e_x (猜測 0.6)
x0_e = [0.4, 0.6];
optimized_e = fminsearch(@(x) cost_function_bounce(x, state_at_bounce_Low, optimized_Cd, dist_L_2nd_actual), x0_e, options);
fprintf('-> 最佳垂直恢復係數 (e_y): %.4f\n', optimized_e(1));
fprintf('-> 最佳水平減速比例 (e_x): %.4f\n', optimized_e(2));
fprintf('\n=== 系統盲點消除，參數鑑別大功告成！ ===\n');

% 副程式區
% 階段一：同時搜尋 Cd 與 仰角 (Theta)
function [best_params, final_state_Low] = fminsearch_Phase1(x0, v0_H, dist_H, v0_L, dist_L, options)
    % 目標函數：高低轉速第一落點誤差的平方和
    cost_func = @(params) (simulate_flight(v0_H, params(1), params(2), 'dist') - dist_H)^2 + ...
                          (simulate_flight(v0_L, params(1), params(2), 'dist') - dist_L)^2;
    best_params = fminsearch(cost_func, x0, options);
    
    % 回傳低轉速落地狀態
    final_state_Low = simulate_flight(v0_L, best_params(1), best_params(2), 'state');
end
% 單次飛行模擬器 (加入仰角三角函數)
function output = simulate_flight(v0, Cd, theta_deg, return_type)
    h0 = 1.21; m = 0.145; rho = 1.225; A = pi * (0.0365)^2; g = 9.81;
    
    % 利用 cosd 與 sind 將初速拆解為 X 與 Y 向量
    vx0 = v0 * cosd(theta_deg);
    vy0 = v0 * sind(theta_deg);
    
    state = [0, h0, vx0, vy0]; 
    dt = 0.01; 
    
    while state(2) > 0
        v = sqrt(state(3)^2 + state(4)^2);
        F_drag = 0.5 * rho * v^2 * Cd * A;
        ax = - (F_drag * (state(3)/v)) / m;
        ay = - g - (F_drag * (state(4)/v)) / m;
        
        state(1) = state(1) + state(3)*dt;
        state(2) = state(2) + state(4)*dt;
        state(3) = state(3) + ax*dt;
        state(4) = state(4) + ay*dt;
    end
    
    if strcmp(return_type, 'dist')
        output = state(1);
    else
        output = state; 
    end
end


% 階段二：計算草地彈跳與第二落點誤差

function SSE = cost_function_bounce(params, initial_state, Cd, dist_2nd_actual)
    e_y = params(1);
    e_x = params(2);
    
    dist_2nd_pred = simulate_second_bounce(initial_state, Cd, e_y, e_x);
    SSE = (dist_2nd_pred - dist_2nd_actual)^2;
end

function dist = simulate_second_bounce(state, Cd, e_y, e_x)
    m = 0.145; rho = 1.225; A = pi * (0.0365)^2; g = 9.81;
    dt = 0.01;
    
    
    state(2) = 0.01; 
    state(4) = -state(4) * e_y; 
    state(3) = state(3) * e_x;  
    
    while state(2) > 0
        v = sqrt(state(3)^2 + state(4)^2);
        F_drag = 0.5 * rho * v^2 * Cd * A;
        ax = - (F_drag * (state(3)/v)) / m;
        ay = - g - (F_drag * (state(4)/v)) / m;
        
        state(1) = state(1) + state(3)*dt;
        state(2) = state(2) + state(4)*dt;
        state(3) = state(3) + ax*dt;
        state(4) = state(4) + ay*dt;
    end
    dist = state(1);
end