# Setting up the Linux environment
* DFT calculations are computational costly, so it is not a good idea to do the calcuation on your desktop/laptop.
* So we usually *submit* our computational procedure as a *job* in the remote environments or supercomputers.
* To use the remote computers, we need to *log-in* from the computer you are using now.
* We use **ssh** communication protocol for this purpose.

* DFTなどの第一原理計算ソフトウェアは計算負荷が大きいので、手元のPCで実行するのは限界がある。そこで、大学の計算機センターやクラウド上のPCにログインして実行する
* これらリモートの環境にログインするには、sshという通信プロトコルを使う。これはlinuxのコマンドになっている
* さらにログインしたリモート環境もlinuxで動くので、結局手元のPCにlinux(っぽい)環境を設定してsshでリモートにログインする

## Mac
* そもそもmacosがlinuxなので特に設定は必要ない
* "アプリケーション" --> "ユーティリティ" --> "ターミナル.app"を開くとsshコマンドが使える

## Windows
* Puttyを使うかWSL2を使う方法がある。ここでは後者を説明する。
1. "スタートボタン"を右クリック --> "Windows Powershell(管理者)"もしくは"Windowsターミナル(管理者)"選ぶ
2. プロンプトで`wsl --install`を入力
  * `wsl --install --online`でインストールできるディストビューションを表示できる
3. 再起動
4. ユーザー名・パスワードを設定

# リモート環境の構築
* sshでログインする環境を設定する。どこかの計算機センターに申し込む必要がある(AWSの無料版とかでもいいのかもしれないがやったことない)
* 東京工業大学のTSUBAMEを使う場合は、別ファイル"tsubame.md"を参照
* 他の大学の計算機センターもやり方は大体同じ
1. ブラウザでポータルサイトにログイン
2. SSH公開鍵の登録を行う部分があるので、そこをクリック
3. ssh-keygenで作成した公開鍵(id_rsa_XXX.pub)を全部コピー
4. 公開鍵を追加する、的な所に貼り付ける
5. ターミナルに戻り、`ssh [your_account_name]@[login-node].ac.jp -i [your_private_keyfile]`でログイン

* ここで、公開鍵は"id_rsa_XXX.pub"、秘密鍵は"id_rsa"といったファイル名になっている
* 5のコマンドは長いので、上記内容を手元のPCの`~/.ssh/config`に登録しておくと次回以降ラクになる
    ```bash
    Host 任意の接続名(asdf)
    HostName ホスト名
    User ユーザー名
    IdentityFile 鍵へのPATH(例えば~/.ssh/asdf.key)
    ```
    * これを作ると`ssh asdf`でログインできる
* 参考: https://qiita.com/soma_sekimoto/items/35845495bc565c38ae9d
