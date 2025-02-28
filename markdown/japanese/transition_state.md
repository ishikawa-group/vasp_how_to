# Transition state (TS, 遷移状態)を求める
* 化学反応において、反応前(reactant、反応物)と反応後(product、生成物)を繋ぐポテンシャルエネルギー曲面上のうち最小エネルギー経路となる経路上でのエネルギー最大値を遷移状態という
* VASPでこのTSを求める
* Nudged elastic band (NEB)という方法がよく用いられる
* Univ. TexasのHenkelman groupのスクリプトがよく用いられる: https://theory.cm.utexas.edu/vtsttools/

* 以下を準備する
1. reactantとproductの構造を求める。各々、構造最適化によって求めること --> POSCAR1, POSCAR2とし同じディレクトリに置く
2. reactantとproductの間の構造を作成する。この中間構造はNEB計算の初期値となるが、適当に両点を補間した構造でよい。ここでは、例として4つの中間構造を作るとする(これをimage数と呼ぶ)
    * 以下を実行: `nebmake.pl POSCAR1 POSCAR2 4`
4. VASPのディレクトリを用意する。4点の場合、`00`, `01` ... `05`が必要。ここで`00`にはreactant、`05`にはproductの構造を置く
5. INCARファイルを次のように作成する
```
IBRION  = 2
POTIM   = 0.1

IMAGES  =  4  # image数に応じて変化させる
SPRING  = -5

ISYM    = 0  # 対称性は切っておいた方がいい
```
5. 計算を実行する
