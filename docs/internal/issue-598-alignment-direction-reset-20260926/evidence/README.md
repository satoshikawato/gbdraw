# 実装証拠の保存方針

Author: `satoshikawato`

各セッションは `S00.md`〜`S04.md` に自分の結果を保存する。未実施の実装・検証を記入済みの成功結果として作らない。

各ファイルに対象commit/base、実行環境、変更したowner/pathと廃止経路、対応する受入条件A01〜A16、fixture、再実行command、結果、既知の制限、次セッションの前提を記す。証拠ファイル自身を含むcommit SHAは自己参照になるため、ファイル内には検証対象runtime commitまたはtree/diff識別子を記し、セッション最終回答で証拠commitとremote SHAを報告する。未変更の証拠は入力・環境・条件が同一の場合だけ再利用する。

画像やログは必要なものだけ保存する。secret、個人のsource data、generated wheel、全量の冗長ログはコミットしない。public screenshotsは実操作と再生成scriptに結び付ける。
