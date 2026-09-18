# ⚾ 智慧棒球守備訓練機 - 系統鑑別與 RPM 控制補償

本專案為大學專題之演算法核心原始碼。主要解決雙輪發球機因空氣阻力與打滑造成的落點誤差，並透過演算法反向抓出硬體未知的隱藏仰角。

## 核心演算法檔案
* `baseball_identification.m`：使用 MATLAB `fminsearch` 進行三變數（空氣阻力 Cd、隱藏仰角 Theta、草地恢復係數 ey/ex）聯合鑑別。
* `Generate_Lookup_Table.m`：導入綜合打滑係數，反向推導 17 公尺實戰守備範圍內之雙落點 RPM 補償資料庫。
