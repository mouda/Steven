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
......winULSA4b2_DC\bin\Release\runHN_FR.95_CR.48_N195_Finsed
......winULSA2i2_MC\bin\Release\runHN_FR.95_CR.48_N195_Finsed
......winULSAkmeans4b2_DC\bin\Release\rrunHN_FR.95_RandMore_CR.48_Finsed
......winULSAkmeans2i_MC\bin\Release\runHN_FR.95CR.48RanMore__Finsed

2.�����k�G
 "DirectAcessPowerMaxFR195_R500_2_MC.txt" �o�� resource consumption �������� (����)�A
 ��Fidelity Ratio=0.95���I�A�Nresource���������Ǥ�ʼg�J Matlab script ���� �ܼ� noClu�C

���U�Ӷi�J�ѤU�� Simulator Path �åB�������� �Ҧ���.bat�ɮ� (�b���S���]�wjobscheduling.bat)
�����᦬�� �A��JData Path������Fig13a.m�Y�i
3.���G��X�G
���X:
ULSA4b2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSA2i2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSAkmeans4b2_All_N195_BW180.0PW0.001_FR0.95.txt
ULSAkmeans2i_All_N195_FR0.95.txt
��Jmatlab plot���|�Y�i�C 

4.�������Ѽ�:
�Цܸ��| �C��simulator path ���W�@�h��Ƨ����A��python��Path�A
�íק�䤤genScriptPy_RunHN.py

5.�Ƶ��G
  K-means�t�C baseline ��Ӧ��ץ�random���覡�A�GK means���G�M�Ϥ��P�A
  �����v�T�۹ﵲ�G
  
  data path ���٦��@�Ӫ����O ULSA2i3 �R�X�Ӫ����G�A�o�O��ӸɤW�����G�A��
  �S���e�X�ӡA�S�I�Oisolate �M join �èS���Ҽ{ correlation �ҥHisolate
  �Mjoin���ĪG�|�ܮt�A���ݭn�i�H�e�e�ݡC
