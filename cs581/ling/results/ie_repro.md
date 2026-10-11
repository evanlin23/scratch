# IE reproduction (RWT screened / unscreened datasets)

## screened: 294 characters (L=259, M=13, P=22), 24 languages; 30 marked "!" (required compatible), 0 marked "*" (removed by screening)
- type L: mean #states 14.41, median 15, max 24
- type M: mean #states 8.77, median 8, max 22
- type P: mean #states 2.32, median 2, max 6
- **Reference tree** compatibility: {'all': (277, 294), 'L': (242, 259), 'M': (13, 13), 'P': (22, 22), 'required(!)': (30, 30), 'extra_steps_total': 19}
  incompatible on reference tree: all1, breast1, float2, head, one, straight, suck2, arm, beard, break1, free, leave1, nine, pour, thousand1, young2, tear

| method | secs | compat (all/L/M/P/!) | RF to ref | Anatolian | Tocharian | Indo-Iranian | Italic | Celtic | Germanic | Balto-Slavic | Baltic | Iranian | Anatolian+Tocharian | Greco-Armenian | Italo-Celtic | Satem core |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 0.0 | 277/242/13/22/30 | 0 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y |
| MP | 4.7 | 280/248/12/20/27 | 10 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | . |
| MC | 0.6 | 280/248/12/20/27 | 10 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | . |
| WMC(!x100) | 0.8 | 278/243/13/22/30 | 4 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y |
| WMC(M,P x5) | 0.8 | 278/243/13/22/30 | 4 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y |
| CapPars(2) | 0.7 | 280/248/12/20/27 | 10 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | . |
| NJ | 0.0 | 276/246/9/21/26 | 12 | Y | Y | Y | Y | Y | Y | Y | Y | Y | . | Y | . | Y |
| ML-Mk | 110.7 | 279/245/12/22/29 | 6 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y |
| ML-binary | 4.8 | 277/243/12/22/29 | 8 | Y | Y | Y | Y | Y | Y | Y | Y | Y | . | Y | Y | Y |

## unscreened: 375 characters (L=336, M=17, P=22), 24 languages; 30 marked "!" (required compatible), 53 marked "*" (removed by screening)
- type L: mean #states 14.27, median 15, max 24
- type M: mean #states 9.94, median 10, max 22
- type P: mean #states 2.32, median 2, max 6
- **Reference tree** compatibility: {'all': (307, 375), 'L': (270, 336), 'M': (15, 17), 'P': (22, 22), 'required(!)': (30, 30), 'extra_steps_total': 84}
  incompatible on reference tree: M9a, M10a, all1, black, blood1, blood2, breast1, fall, father, fire, float2, fog1, fog2, head, hear1, hear2, heavy2, here, hold, horn1, horn2, husband, I, me, if, lie, neck2, one, sleep, snake1, snake2, snow1, snow2, straight, suck2, swim, thee, tooth, we, where, ye, arm, arrow, beard, break1, bull1, bull2, duck, free, grain2, honey, house1, house2, lamb, leave1, nine, ox, pig2, pour, put1, stay1, stay2, sweat, thousand1, weave, wolf, young2, tear

| method | secs | compat (all/L/M/P/!) | RF to ref | Anatolian | Tocharian | Indo-Iranian | Italic | Celtic | Germanic | Balto-Slavic | Baltic | Iranian | Anatolian+Tocharian | Greco-Armenian | Italo-Celtic | Satem core |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| reference | 0.0 | 307/270/15/22/30 | 0 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y |
| MP | 0.9 | 317/282/15/20/27 | 10 | Y | Y | Y | Y | Y | Y | Y | Y | Y | . | Y | Y | . |
| MC | 1.1 | 319/284/15/20/27 | 8 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | . |
| WMC(!x100) | 1.1 | 308/271/15/22/30 | 2 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y |
| WMC(M,P x5) | 1.0 | 317/281/16/20/27 | 8 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | . |
| CapPars(2) | 1.0 | 319/284/15/20/27 | 8 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | . |
| NJ | 0.0 | 313/277/16/20/27 | 10 | Y | Y | Y | Y | Y | Y | Y | Y | Y | . | Y | Y | . |
| ML-Mk | 175.9 | 311/275/14/22/29 | 6 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y |
| ML-binary | 6.5 | 310/274/14/22/29 | 6 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y |

