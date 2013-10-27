Fig14a
Author:Steven

1.Background：
...Simulator: 
......DirectAcessPowerMax_MC (與Base Station直接傳輸)
......ULSA4b2 (Two-Tier Data-Centric Clustering) 
......ULSA3i2 (Two-Tier Data-Centric Clustering (with no isolate and join))
...Simulator Path:
......DirectAcessPowerMax_MC\bin\Release
......winULSA4b2_DC\bin\Release\runHN_FR.95_CR.48_N195_Finsed
......winULSA3i2_DC\bin\Release\runHN_FR.95_CR.48_N195

2.執行方法：
 "DirectAcessPowerMaxFR195_R500_2_MC.txt" 得到 resource consumption 的比較基準 (分母)，
 拿Fidelity Ratio=0.95的點，將resource視為比較基準手動寫入 Matlab script 中的 變數 noClu。

接下來進入剩下的 Simulator Path 執行 在jobscheduling.bat 結束後收割 ，放入Data Path中執行Fig14a.m即可
3.結果輸出：
拿出:
ULSA4b2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSA3i2_All_N195_BW180.0PW0.001_FR0.95.txt
放入matlab plot路徑即可。 

4.更改模擬參數:
請至路徑 每個simulator path 的上一層資料夾中，找python的Path，
並修改其中genScriptPy_RunHN.py

5.備註：
ULSA4b2_All_N195_BW180.0PW0.001_FR0.95.txt 是與 Fig 13a 用同一份答案 
