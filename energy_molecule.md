# エネルギー計算
* まずCO分子のエネルギーを計算してみる
* ここで言うエネルギーとは、Schorodinger方程式から得られる"電子エネルギー + 核間反発エネルギー"のこと
* 必要なファイルは`INCAR`, `POSCAR`, `KPOINTS`, `POTCAR`の4つなのでこれを作成する
* これらを新規ディレクトリに入れて、そのディレクトリ内で実行することにする
* (もしVASPの実行ファイルが手元のPCにあるのであれば)上記のファイルがあるディレクトリで`vasp`とすれば、計算が実行される(`vasp`のパスが通っている場合)。ただし、ほとんどの場合はjobとして実行するので、このように実行することはあまりない
* 本例題はvaspwiki https://www.vasp.at/wiki/index.php/CO を参考にした

## INCAR
* VASPの計算条件を指示するファイルで、キーワードと変数のペアを入力する
* 変数は、integer(整数)、float(浮動小数点実数)、logical(.TRUE. or .FALSE.)などがある
* キーワードの意味については`incar.md`で詳しく解説するので、今はとりあえず下のようなファイルとする
```
SYSTEM = CO molecule in a box
ISMEAR = 0  ! Gaussian smearing
NSW    = 5  ! 5 ionic steps
IBRION = 2  ! use the conjugate gradient algorithm
```

## POSCAR
* 原子の位置を与えるファイル
* 単位セル(unit cell)の大きさ、元素名とその個数、各元素のXYZ座標という順のセクションとなっている
```
CO molecule in a box
 1.0          ! universal scaling parameters
 8.0 0.0 0.0  ! lattice vector  a(1)
 0.0 8.0 0.0  ! lattice vector  a(2)
 0.0 0.0 8.0  ! lattice vector  a(3)
1 1           ! number of atoms for each species
cart          ! positions in cartesian coordinates
 0 0 0        ! first atom
 0 0 1.12     ! second atom
 ```

## KPOINTS
* 逆格子空間を数値積分で求める際の、k点の個数と位置を記述する
```
Gamma-point only
 0
Monkhorst Pack
 1 1 1
 0 0 0
 ```
 
## POTCAR
* 擬ポテンシャル(pseudopotential)の情報を書き込む
* VASPでは、特にprojector augumented-wave (PAW)という呼び方をする
* VASPのインストールディレクトリにあるファイルから以下のコマンドで生成する
`cat  .../O/POTCAR  .../C POTCAR  > POTCAR`
 * 上記が準備できたら、`vasp_std`を実行

# ジョブの実行
* 以下ではTSUBAMEを想定する
* 計算機センター等ではjob queueing systemというのをよく使う。これは、空いている計算機にジョブを実行させるためのシステムである
* ログインノードでジョブをsubmitして、計算ノードで計算を実行する。したがって計算ノードにログインする必要がない
* 何をして欲しいかをジョブスクリプトと呼ばれるファイル(.sh形式)に記述して、それを計算ノードに実行してもらう
* したがってジョブスクリプトファイルを書く必要がある
* 本格的な計算をするには、計算資源を購入して指定する必要がある。以下では、「お試しモード」というのがあるのでそちらを利用する
* お試しモードは課金なしだが制限時間10分。あくまで計算がうまく流れるかの確認用なので注意。

## ジョブスクリプト
* MPIと呼ばれる並列計算を可能にする機能を利用する
* https://helpdesk.t3.gsic.titech.ac.jp/manuals/handbook.ja/jobs/ のMPI並列、intel-mpi利用をベースにVASP実行用に編集する
* 以下のようなファイルを`run.sh`を作成する。詳細は意味については、今はわからなくても大丈夫

```bash
#!/bin/sh
#$ -cwd
# 資源タイプF 1ノードを使用
#$ -l f_node=1
#$ -l h_rt=0:10:00
#$ -N flatmpi

# VASPのパスを指定
PRG=/home/your_vasp_dir/bin/vasp_std

. /etc/profile.d/modules.sh
module load cuda
module load intel

# Intel MPI環境の読込
module load intel-mpi

# ノードあたり8プロセスMPI全32 プロセスを使用
mpiexec.hydra -ppn 8 -n 8 ${PRG} >& vasp.out
```

## ジョブの投入
* qsub: `qsub -g [TSUBAME_group] run.sh`
  * TSUBAME_group: 自分のTSUBAMEグループ名。ただしお試しモードでは必要ないので書かないことにする
* 上記までがVASPのジョブ投入になる。うまく実行できていると、アウトプットファイル(OUTCARなど)が出てくる

## 以下は必要に応じて使う
### ジョブの状態確認
* qstatコマンドでジョブの状態を確認できる: `qstat [option]`

### ジョブの中止
* qdel: `qdel JOB_ID`
