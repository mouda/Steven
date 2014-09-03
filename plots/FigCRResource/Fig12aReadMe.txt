Fig12a
Author:Steven

1.Background：
...Simulator: 
              DirectAcessPowerMax_MC (與Base Station直接傳輸)
              ULSA4b2 (Two-Tier Data-Centric Clustering) 
              ULSA2i2 (Two-Tier Machine-Centric Clustering)
              winULSAkmeans4b2_DC (Two-Tier Data-Centric K-Means) 
              winULSAkmeans2i_MC (Two-Tier Machine-Centric K-Means)
...Simulator Path:
......DirectAcessPowerMax_MC\bin\Release
......winULSA4b2_DC\bin\Release\runCR_FR.95_Finsed
......winULSA2i2_MC\bin\Release\runCR_FR.95_Finsed
......winULSAkmeans4b2_DC\bin\Release\runCR_FR.95_Finsed
......winULSAkmeans2i_MC\bin\Release\runCR_FR.95_Finsed

2.執行方法：
先執行 DirectAcessPowerMax_MC路徑中的RunByCRDirect.bat，即可從 
"DirectAcessPowerMaxFR195_R500_2_MC.txt" 得到 resource consumption 的比較基準 (分母)，
將比較基準寫入 (一個 Array) 手動寫入 Matlab script 中的 變數 noClu。

接下來進入剩下的 Simulator Path 並且直接執行 所有的.bat檔案 (在此沒有設定jobscheduling.bat)
結束後收割 ，放入Data Path中執行Fig12a.m即可
3.結果輸出：
拿出:
ULSA4b2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSA2i2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSAkmeans4b2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSAkmeans2i_All_N195_FR0.95.txt
放入matlab plot路徑即可。 

4.更改模擬參數:
請至路徑 每個simulator path 的上一層資料夾中，找python的Path，
並修改其中genScriptPy_RunCR.py

