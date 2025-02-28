# 東京工業大学の計算機センター(TSUBAME4.0)を利用するためのメモ
* 以下、東工大の職員・学生を想定

## アカウントの作成
1. 東工大ポータルにログイン、"TSUBAMEポータル"があるのでそこをクリック
2. 必要事項を入力
3. 利用の承認があると、メールが来るのでそこのURLをクリック
4. アカウント名を確認、SSH公開鍵の登録を行う
    4-1. mac/windows-wslのターミナルを立ち上げて、`ssh-keygen`を入力
    4-2. 公開鍵のファイル名を聞かれるので、適宜入力
    4-3. 公開鍵ファイル`id_rsa_XXX_.pub`の内容をコピー
    4-4. TSUBAMEポータルに戻り、"SSH公開鍵追加"をクリック、公開鍵コード入力に上記の内容を貼り付ける
    4-5. "追加"をクリック
5. ターミナルに戻り、`ssh [your_account_name]@login.t4.gsic.titech.ac.jp -i [your_private_keyfile]`でログイン
6. ssh configを設定

## リソースの獲得
* 計算機リソースを得るためにはグループを作る必要がある

#### グループの作成
1. TSUBAMEポータルにログイン
2. "グループ作成"をクリック、必要事項を入力するとグループができる
3. 作成後、"所属グループ管理" --> "詳細表示"でグループメンバーを増やしたり減らしたりできる

#### 課金について
##### 支払いコード登録
1. "支払いコード"が必要になるのでそれを作る
2. "支払いコード管理"をクリック --> "新規支払いコード申請"
3. 必要事項を入力。予算コードや予算名称が必要になる。
4. 承認されるのを待つ
5. 予算コードが発行されたら、TSUBAMEポータルの"所属グループ管理"から"TSUBAMEポイント"の部分があるので、そこで課金情報を入力する

##### ポイント購入
1. TSUBAMEポータルにログイン
2. グループ管理
3. 所属グループを選択
4. 詳細
5. ポイント購入

#### TSUBAMEポイント確認
* ターミナルで`t4-user-info group point`

## submit
* `qsub -g [group_name] script.sh`
* グループ名を指定しないとお試し実行になる

## ジョブ確認
qstat

## Gaussian
* `module load gaussian`
* サンプルsh
```bash
#!/bin/bash
#$ -cwd
#$ -l f_node=1
#$ -l h_rt=00:10:00
#$ -V

. /etc/profile.d/modules.sh
module load gaussian16

g16 glycine.com
```

### Gaussview
* `qrsh -g tga-ishikawalab -l s_core=1 -l h_rt=1:00:00`でinteractive jobを立ち上げ
* `module load gaussview`して`gview`で立ち上がる
* X経由なのでとても遅い
