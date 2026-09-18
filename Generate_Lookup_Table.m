
% 棒球守備訓練機 - RPM 補償對照表 (Lookup Table) 雙落點終極生成器
% 包含：系統鑑別參數導入、雙轉速綜合打滑補償、第一/第二落點同步追蹤

clear; clc;

%% 1. 匯入已鑑別完成的系統參數 (Ground Truth 演算結果)
Cd = 0.3976;       % 鑑別出的空氣阻力係數
theta_deg = 5.13;  % 鑑別出的機台隱藏仰角 (度)
e_y = 0.4370;      % 垂直恢復係數
e_x = 0.9683;      % 水平減速比例

%% 2. 機器硬體設定與綜合打滑係數計算
r_roller = 0.147;  % 滾輪半徑 14.7 cm (0.147 m)

% --- 高轉速組實驗數據 ---
rpm_H = 1400;         
v0_actual_H = 18.67;  
v0_ideal_H = rpm_H * 2 * pi / 60 * r_roller; 
slip_H = v0_actual_H / v0_ideal_H; 

% --- 低轉速組實驗數據 ---
rpm_L = 1000;         
v0_actual_L = 13.14;  
v0_ideal_L = rpm_L * 2 * pi / 60 * r_roller; 
slip_L = v0_actual_L / v0_ideal_L; 

% --- 綜合打滑係數 (取高低轉速之平均) ---
slip_ratio = (slip_H + slip_L) / 2;

fprintf('=== 系統參數與硬體特性 ===\n');
fprintf('空氣阻力 Cd: %.4f, 仰角: %.2f 度\n', Cd, theta_deg);
fprintf('草地彈性 e_y: %.4f, 摩擦力 e_x: %.4f\n', e_y, e_x);
fprintf('高轉速打滑係數: %.4f, 低轉速打滑係數: %.4f\n', slip_H, slip_L);
fprintf('=> 系統綜合打滑係數 (Slip Ratio) 鎖定為: %.4f\n\n', slip_ratio);

%% 3. 設定要建立對照表的「目標守備距離」

target_distances = 5:1:17; 
target_type = '2nd'; % 鎖定以 '2nd' (第二彈跳) 為主要目標距離

%% 4. 反向推導 RPM (建立雙落點 Lookup Table)
fprintf('=== 棒球發球機 RPM 補償對照表 (雙落點追蹤) ===\n');
fprintf('第二落點(目標)\t第一落點(預估)\t所需初速(m/s)\t馬達設定(RPM)\n');
fprintf('-------------------------------------------------------------\n');

% 準備一個陣列來儲存對照表 [目標距離, 第一落點, 初速, RPM]
lookup_table = zeros(length(target_distances), 4);
options = optimset('Display', 'off');

for i = 1:length(target_distances)
    target_dist = target_distances(i);
    
    % 1. 使用最佳化演算法，反向找尋「為了達到第二落點，需要多快的初速 v0」
    v0_required = fminsearch(@(v0) cost_func_find_v0(v0, target_dist, target_type, Cd, theta_deg, e_y, e_x), 15, options);
    
    % 2. 初速確定後，順便把這個初速帶回模擬器，算出「第一落點」在哪裡
    dist_1st_pred = simulate_flight_to_bounce(v0_required, Cd, theta_deg);
    
    % 3. 考慮打滑係數，把初速反推回馬達必須輸出的「理想速度」
    v0_ideal_required = v0_required / slip_ratio;
    
    % 4. 將理想速度轉換回馬達面板該設定的 RPM
    rpm_required = v0_ideal_required / r_roller * 60 / (2*pi);
    
    % 記錄並印出結果
    lookup_table(i, :) = [target_dist, dist_1st_pred, v0_required, rpm_required];
    fprintf('%2d m\t\t\t%5.2f m\t\t\t%5.2f\t\t\t%4.0f\n', target_dist, dist_1st_pred, v0_required, round(rpm_required));
end

%% 5. 繪製補償對照表曲線圖 (雙軌跡視覺化)
figure;
% 畫出第二落點 (主要目標) 曲線
plot(lookup_table(:,1), lookup_table(:,4), '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', '目標落點 (第二彈跳)');
hold on;
% 畫出第一落點 (參考預估) 曲線
plot(lookup_table(:,2), lookup_table(:,4), '-s', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'Color', 'r', 'DisplayName', '預估落點 (第一彈跳)');

xlabel('落點距離 (公尺)');
ylabel('發球機馬達設定轉速 (RPM)');
title(sprintf('棒球發球機 RPM 補償對照表 (機台仰角: %.2f度)', theta_deg));
legend('Location', 'southeast');
grid on;
hold off;



% 副程式區：包含飛行與彈跳的物理模擬核心

function SSE = cost_func_find_v0(v0, target_dist, target_type, Cd, theta_deg, e_y, e_x)
    % 根據輸入的速度 v0，模擬其落點，並與目標距離計算誤差
    if strcmp(target_type, '1st')
        dist_pred = simulate_flight_to_bounce(v0, Cd, theta_deg);
    else
        state_at_bounce = simulate_flight_to_bounce_state(v0, Cd, theta_deg);
        dist_pred = simulate_second_bounce(state_at_bounce, Cd, e_y, e_x);
    end
    SSE = (dist_pred - target_dist)^2;
end

function dist = simulate_flight_to_bounce(v0, Cd, theta_deg)
    h0 = 1.21; m = 0.145; rho = 1.225; A = pi * (0.0365)^2; g = 9.81;
    state = [0, h0, v0*cosd(theta_deg), v0*sind(theta_deg)]; dt = 0.01; 
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

function state_out = simulate_flight_to_bounce_state(v0, Cd, theta_deg)
    h0 = 1.21; m = 0.145; rho = 1.225; A = pi * (0.0365)^2; g = 9.81;
    state = [0, h0, v0*cosd(theta_deg), v0*sind(theta_deg)]; dt = 0.01; 
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
    state_out = state;
end

function dist = simulate_second_bounce(state, Cd, e_y, e_x)
    m = 0.145; rho = 1.225; A = pi * (0.0365)^2; g = 9.81; dt = 0.01;
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