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
......winULSA4b2_DC\bin\Release\runHN_FR.95_CR.48_N195_Finsed
......winULSA2i2_MC\bin\Release\runHN_FR.95_CR.48_N195_Finsed
......winULSAkmeans4b2_DC\bin\Release\rrunHN_FR.95_RandMore_CR.48_Finsed
......winULSAkmeans2i_MC\bin\Release\runHN_FR.95CR.48RanMore__Finsed

2.執行方法：
 "DirectAcessPowerMaxFR195_R500_2_MC.txt" 得到 resource consumption 的比較基準 (分母)，
 拿Fidelity Ratio=0.95的點，將resource視為比較基準手動寫入 Matlab script 中的 變數 noClu。

接下來進入剩下的 Simulator Path 並且直接執行 所有的.bat檔案 (在此沒有設定jobscheduling.bat)
結束後收割 ，放入Data Path中執行Fig13a.m即可
3.結果輸出：
拿出:
ULSA4b2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSA2i2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSAkmeans4b2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSAkmeans2i_All_N195_FR0.95.txt
放入matlab plot路徑即可。 

4.更改模擬參數:
請至路徑 每個simulator path 的上一層資料夾中，找python的Path，
並修改其中genScriptPy_RunHN.py

5.備註：
  K-means系列 baseline 後來有修正random的方式，故K means結果和圖不同，
  但不影響相對結果
  
  data path 中還有一個版本是 ULSA2i3 吐出來的結果，這是後來補上的結果，並
  沒有畫出來，特點是isolate 和 join 並沒有考慮 correlation 所以isolate
  和join的效果會變差，有需要可以畫畫看。
