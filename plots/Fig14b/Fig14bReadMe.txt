Fig14b
Author:Steven

1.Background：
...Simulator: 
              ULSA4b2 (Two-Tier Data-Centric Clustering) 
...Simulator Path:
......winULSA4b2_DC\bin\Release\runHN_FR.95_CR.48_N50_Finsed
......winULSA4b2_DC\bin\Release\runHN_FR.95_CR.48_N150_Finsed
......winULSA4b2_DC\bin\Release\runHN_FR.95_CR.48_N195_Finsed

2.執行方法：
進入Simulator Path 執行在jobscheduling.bat，在jobscheduling.bat 
結束後收割 ，放入Data Path中執行Fig14b.m即可
3.結果輸出：
拿出:
ULSA4b2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSA4b2_All_N150_BW180.0PW0.001_FR0.95.txt
ULSA4b2_All_N50_BW180.0PW0.001_FR0.95.txt
放入matlab plot路徑即可。 

4.更改模擬參數:
請至路徑 每個simulator path 的上一層資料夾中，找python的Path，
並修改其中genScriptPy_RunHN.py

