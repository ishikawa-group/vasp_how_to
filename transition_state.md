# Transition state (TS, 遷移状態)を求める
* 化学反応において、反応前(reactant、反応物)と反応後(product、生成物)を繋ぐポテンシャルエネルギー曲面上のうち最小エネルギー経路となる経路上でのエネルギー最大値を遷移状態という
* VASPでこのTSを求める
* Nudged elastic band (NEB)という方法がよく用いられる
* Univ. TexasのHenkelman groupのスクリプトがよく用いられる

* 以下を準備する
1. reactantとproductの構造を求める。各々、構造最適化によって求めること
2. reactantとproductの間の構造を作成する。この中間構造はNEB計算の初期値となるが、適当に両点を補間した構造でよい。ここでは、例として4つの中間構造を作るとする
3. VASPのディレクトリを用意する。4点の場合、`00`, `01` ... `05`が必要。ここで`00`にはreactant、`05`にはproductの構造を置く
4. INCARファイルを次のように作成する
```
 ISTART = 1
 ICHARG = 1
 INIWAV = 1
 ISPIN  = 1
 ISYM   = 0      !<=== おそらく対称性の制限を外した方が良い

 ENCUT  =  400.000
 LREAL  =  Auto
 EDIFF  =  1.0e-04
 ALGO   =  Fast
 PREC   =  Normal

 NELMIN =  4
 NELM   =  100
 NBANDS =   98

 ISMEAR =   2
 SIGMA  =   0.2

 IBRION =    1    !<=== RMM-DIIS (1) or quick-min (3) が推奨されている
 NSW    =   100
 ISIF   =    2

 NCORE  =    8

 # NEB OPTIONS
 IMAGES = 7       !<=== 始点と終点を除くイメージ数
 SPRING = -5.0    !<=== バネ定数
 LCLIMB = .TRUE.  !<=== climbing image NEBを行うかどうか
```
5. 計算を実行する
