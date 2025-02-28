# 計算結果の表示・確認・可視化
* 計算が終了したら、その結果を可視化したい。ここではその方法について説明する

## 構造

## 状態密度
* orbital?が各エネルギー幅にどのくらい詰まっているかを表したもの
* density of state (DOS)という
* DOSの値が大きいほど、そのエネルギー領域に密に状態が詰まっていることを意味する
* k点の値がある程度大きくないと正しいDOSにならないことに注意
* INCARのsmearingパラメーター(ISMEAR, SIGMA)もDOSの精度に関係する

* DOSCARというファイルが出力される
* spin polarizedとnon-polarizedで出力が違うので注意
* non-polarizedの場合、各項目の意味は
  * `energy DOS integrated-DOS`
* 原子ごとに別々にしたDOSはlocalized DOS (LDOS)という
* これに対して、s, p, dなどの角運動量に分解したものはprojected DOS (PDOS)という
* PDOSを得るためにはINCARでLORBIT=11としておかなければならないので注意

## 電荷密度
* VESTAとうプログラムを使うのが良い
* CHGCARをダウンロードし、VESTAからロードすれば見れる
* VESTAでの電子密度の単位は$1/Bohr^3$となっている
