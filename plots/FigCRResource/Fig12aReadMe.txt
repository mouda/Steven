Fig12a
Author:Steven

1.Background�G
...Simulator: 
              DirectAcessPowerMax_MC (�PBase Station�����ǿ�)
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

2.�����k�G
������ DirectAcessPowerMax_MC���|����RunByCRDirect.bat�A�Y�i�q 
"DirectAcessPowerMaxFR195_R500_2_MC.txt" �o�� resource consumption �������� (����)�A
�N�����Ǽg�J (�@�� Array) ��ʼg�J Matlab script ���� �ܼ� noClu�C

���U�Ӷi�J�ѤU�� Simulator Path �åB�������� �Ҧ���.bat�ɮ� (�b���S���]�wjobscheduling.bat)
�����᦬�� �A��JData Path������Fig12a.m�Y�i
3.���G��X�G
���X:
ULSA4b2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSA2i2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSAkmeans4b2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSAkmeans2i_All_N195_FR0.95.txt
��Jmatlab plot���|�Y�i�C 

4.�������Ѽ�:
�Цܸ��| �C��simulator path ���W�@�h��Ƨ����A��python��Path�A
�íק�䤤genScriptPy_RunCR.py

