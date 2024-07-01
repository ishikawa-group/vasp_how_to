# VASPのコンパイル
* Vienna ab initio simulation package (VASP)をリモート環境でコンパイルする
* コンパイルとは、ソースコードからバイナリーファイル(実行可能ファイル)を生成すること
* 実行可能ファイルが用意されていて、それを使う場合はこの項目は飛ばしてよい

## ソースコードの入手
* VASPのソースコードを入手する
1. VASP社に料金を支払ってVASPのポータルサイトのアカウントを取得
2. ポータルサイトからvaspとpotentialの.tar.gzファイルをダウンロードしてくる

## コンパイル
1. リモート環境にvaspとpotentialのtar.gzをコピーする: `scp vasp.5.4.4.tar.gz remote:/home/your_vasp_dir`
  * remoteはリモートホスト名。`~/.ssh/config`のHostに書いてあるもの。
2. ファイルが送れたらリモート環境に移動
3. vaspを解凍: `tar zxvf vasp.5.4.4.tar.gz`
4. potentialを解凍
  ```
  mkdir -p poentials/potpaw_PBE.54
  cd potentials/potpaw_PBE.54
  mv ../../potpaw_PBE.54.tar.gz ./
  tar zxvf potpaw_PBE.54.tar.gz
  ```
4. `cd vasp.5.4.4`
5. intelコンパイラを利用可能にする: `module load intel; module load intel-mpi`
6. makefileをコピー
    * Intel one apiコンパイラ: `cp arch/makefile.include.oneapi  ./makefile.include`
    * Intelコンパイラ(古め): `cp arch/makefile.include.linux_intel  ./makefile.include`
7. コンパイル: `make all`
8. `bin`ディレクトリに`vasp_std`ができる。これがVASPの実行ファイルになる

### VTSTのコンパイル
* University of Texas, Henkelman groupで開発されているコード
* VASPに追加する形でCINEB, dimer法などの機能を加えることができる
* ソースコードで配布されているので、VASPをコンパイルし直す必要がある

1. vtstcode.tar.gzを入手する: https://theory.cm.utexas.edu/vtsttools/download.html
2. vtstcode6.4のファイルをsrc以下にコピーする
  * vasp6.4.2ディレクトリで`cp ../vtstcode.6.4.2/* ./src`
3. `src/main.F`の三箇所を以下のように書き換える
  a.
    ```fortran
    CALL CHAIN_FORCE(T_INFO%NIONS,DYN%POSION,TOTEN,TIFOR, &
        LATT_CUR%A,LATT_CUR%B,IO%IU6)
    ```
    ↓
    ```fortran
    CALL CHAIN_FORCE(T_INFO%NIONS,DYN%POSION,TOTEN,TIFOR, &
         TSIF,LATT_CUR%A,LATT_CUR%B,IO%IU6)
    ```
    b.
    ```fortran
    IF (LCHAIN) CALL chain_init( T_INFO, IO)
    ```
    ↓
    ```fortran
    CALL chain_init( T_INFO, IO)
    ```
    c.
    ```fortran
    IF ( DYN%ISIF >= 3 ) THEN
        CALL CHAIN_STRESS( TSIF )
    END IF
    ```
    ↓
    ```fortran
    ! IF ( DYN%ISIF >= 3 ) THEN
    !     CALL CHAIN_STRESS( TSIF )
    ! END IF
    ```
4. `src/.objects`の**chain.oの前**に以下を追加する
    ```
    bfgs.o dynmat.o instanton.o lbfgs.o sd.o cg.o dimer.o bbm.o \
    fire.o lanczos.o neb.o qm.o \
    pyamff_fortran/*.o ml_pyamff.o \
    opt.o \
    ```
5. `src/makefile`の以下を変更する
    ```
    LIB=lib parser
    ```
    ↓
    ```
    LIB= lib parser pyamff_fortran
    ```
6. `make all`

* vaspの元々のソースコードの`main.F`にある`subroutine chain_stress`をvtstの`chain.F`に貼り付けてコンパイルすることもできるので、こちらのほうが良いかもしれないが`ISIF>3`でMD計算をしなければ問題ないと思う
* 参考: https://theory.cm.utexas.edu/vtsttools/installation.html
