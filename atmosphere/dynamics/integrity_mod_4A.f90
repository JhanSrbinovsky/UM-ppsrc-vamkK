! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! tool to test integrity of fields to trace memory errors and to
! compare runs to bit-level
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Dynamics

MODULE integrity_mod

USE chsunits_mod, ONLY : nunits
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

INTEGER, PARAMETER       :: registry_length = 200
CHARACTER (len=5),  SAVE :: key(registry_length)
CHARACTER (len=65), SAVE :: hash(registry_length+1)

INTEGER, SAVE            :: last_entry=0

INTEGER, SAVE            :: av_out=0

LOGICAL, SAVE            :: integrity_test = .FALSE.
LOGICAL, SAVE            :: integrity_test_ghash = .FALSE.
LOGICAL, SAVE            :: integrity_file       = .FALSE.
INTEGER, SAVE            :: global_unt = 0

CONTAINS

SUBROUTINE check_hash_m(field ,field_size ,key_to_check,&
                       field1,field_size1,key_to_check1,&
                       field2,field_size2,key_to_check2,&
                       field3,field_size3,key_to_check3,&
                       field4,field_size4,key_to_check4,&
                       field5,field_size5,key_to_check5,&
                       field6,field_size6,key_to_check6,&
                       field7,field_size7,key_to_check7,&
                       field8,field_size8,key_to_check8,&
                       field9,field_size9,key_to_check9,&
                       field10 ,field_size10 ,key_to_check10,&
                       field11 ,field_size11 ,key_to_check11,&
                       field12 ,field_size12 ,key_to_check12,&
                       field13 ,field_size13 ,key_to_check13,&
                       field14 ,field_size14 ,key_to_check14,&
                       field15 ,field_size15 ,key_to_check15,&
                       field16 ,field_size16 ,key_to_check16,&
                       field17 ,field_size17 ,key_to_check17,&
                       field18 ,field_size18 ,key_to_check18,&
                       field19 ,field_size19 ,key_to_check19,&
                       field20 ,field_size20 ,key_to_check20,&
                       field21 ,field_size21 ,key_to_check21,&
                       field22 ,field_size22 ,key_to_check22,&
                       field23 ,field_size23 ,key_to_check23,&
                       field24 ,field_size24 ,key_to_check24,&
                       field25 ,field_size25 ,key_to_check25,&
                       field26 ,field_size26 ,key_to_check26,&
                       field27 ,field_size27 ,key_to_check27,&
                       field28 ,field_size28 ,key_to_check28,&
                       field29 ,field_size29 ,key_to_check29,&
                       field30 ,field_size30 ,key_to_check30,&
                       field31 ,field_size31 ,key_to_check31,&
                       field32 ,field_size32 ,key_to_check32,&
                       field33 ,field_size33 ,key_to_check33,&
                       field34 ,field_size34 ,key_to_check34,&
                       field35 ,field_size35 ,key_to_check35,&
                       field36 ,field_size36 ,key_to_check36,&
                       field37 ,field_size37 ,key_to_check37,&
                       field38 ,field_size38 ,key_to_check38,&
                       field39 ,field_size39 ,key_to_check39,&
                       field40 ,field_size40 ,key_to_check40,&
                       field41 ,field_size41 ,key_to_check41,&
                       field42 ,field_size42 ,key_to_check42,&
                       field43 ,field_size43 ,key_to_check43,&
                       field44 ,field_size44 ,key_to_check44,&
                       field45 ,field_size45 ,key_to_check45,&
                       field46 ,field_size46 ,key_to_check46,&
                       field47 ,field_size47 ,key_to_check47,&
                       field48 ,field_size48 ,key_to_check48,&
                       field49 ,field_size49 ,key_to_check49,&
                       field50 ,field_size50 ,key_to_check50,&
                       field51 ,field_size51 ,key_to_check51,&
                       field52 ,field_size52 ,key_to_check52,&
                       field53 ,field_size53 ,key_to_check53,&
                       field54 ,field_size54 ,key_to_check54,&
                       field55 ,field_size55 ,key_to_check55,&
                       field56 ,field_size56 ,key_to_check56,&
                       field57 ,field_size57 ,key_to_check57,&
                       field58 ,field_size58 ,key_to_check58,&
                       field59 ,field_size59 ,key_to_check59,&
                       field60 ,field_size60 ,key_to_check60,&
                       field61 ,field_size61 ,key_to_check61,&
                       field62 ,field_size62 ,key_to_check62,&
                       field63 ,field_size63 ,key_to_check63,&
                       field64 ,field_size64 ,key_to_check64,&
                       field65 ,field_size65 ,key_to_check65,&
                       field66 ,field_size66 ,key_to_check66,&
                       field67 ,field_size67 ,key_to_check67,&
                       field68 ,field_size68 ,key_to_check68,&
                       field69 ,field_size69 ,key_to_check69,&
                       field70 ,field_size70 ,key_to_check70,&
                       field71 ,field_size71 ,key_to_check71,&
                       field72 ,field_size72 ,key_to_check72,&
                       field73 ,field_size73 ,key_to_check73,&
                       field74 ,field_size74 ,key_to_check74,&
                       field75 ,field_size75 ,key_to_check75,&
                       field76 ,field_size76 ,key_to_check76,&
                       field77 ,field_size77 ,key_to_check77,&
                       field78 ,field_size78 ,key_to_check78,&
                       field79 ,field_size79 ,key_to_check79,&
                       field80 ,field_size80 ,key_to_check80,&
                       field81 ,field_size81 ,key_to_check81,&
                       field82 ,field_size82 ,key_to_check82,&
                       field83 ,field_size83 ,key_to_check83,&
                       field84 ,field_size84 ,key_to_check84,&
                       field85 ,field_size85 ,key_to_check85,&
                       field86 ,field_size86 ,key_to_check86,&
                       field87 ,field_size87 ,key_to_check87,&
                       field88 ,field_size88 ,key_to_check88,&
                       field89 ,field_size89 ,key_to_check89,&
                       field90 ,field_size90 ,key_to_check90,&
                       field91 ,field_size91 ,key_to_check91,&
                       field92 ,field_size92 ,key_to_check92,&
                       field93 ,field_size93 ,key_to_check93,&
                       field94 ,field_size94 ,key_to_check94,&
                       field95 ,field_size95 ,key_to_check95,&
                       field96 ,field_size96 ,key_to_check96,&
                       field97 ,field_size97 ,key_to_check97,&
                       field98 ,field_size98 ,key_to_check98,&
                       field99 ,field_size99 ,key_to_check99) 

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (len=5)  key_to_check
INTEGER            field_size
REAL               field(*)

CHARACTER (len=5), OPTIONAL ::  key_to_check1
INTEGER          , OPTIONAL ::  field_size1
REAL             , OPTIONAL ::  field1(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check2
INTEGER          , OPTIONAL ::  field_size2
REAL             , OPTIONAL ::  field2(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check3
INTEGER          , OPTIONAL ::  field_size3
REAL             , OPTIONAL ::  field3(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check4
INTEGER          , OPTIONAL ::  field_size4
REAL             , OPTIONAL ::  field4(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check5
INTEGER          , OPTIONAL ::  field_size5
REAL             , OPTIONAL ::  field5(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check6
INTEGER          , OPTIONAL ::  field_size6
REAL             , OPTIONAL ::  field6(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check7
INTEGER          , OPTIONAL ::  field_size7
REAL             , OPTIONAL ::  field7(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check8
INTEGER          , OPTIONAL ::  field_size8
REAL             , OPTIONAL ::  field8(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check9
INTEGER          , OPTIONAL ::  field_size9
REAL             , OPTIONAL ::  field9(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check10
INTEGER          , OPTIONAL ::  field_size10
REAL             , OPTIONAL ::  field10(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check11
INTEGER          , OPTIONAL ::  field_size11
REAL             , OPTIONAL ::  field11(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check12
INTEGER          , OPTIONAL ::  field_size12
REAL             , OPTIONAL ::  field12(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check13
INTEGER          , OPTIONAL ::  field_size13
REAL             , OPTIONAL ::  field13(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check14
INTEGER          , OPTIONAL ::  field_size14
REAL             , OPTIONAL ::  field14(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check15
INTEGER          , OPTIONAL ::  field_size15
REAL             , OPTIONAL ::  field15(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check16
INTEGER          , OPTIONAL ::  field_size16
REAL             , OPTIONAL ::  field16(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check17
INTEGER          , OPTIONAL ::  field_size17
REAL             , OPTIONAL ::  field17(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check18
INTEGER          , OPTIONAL ::  field_size18
REAL             , OPTIONAL ::  field18(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check19
INTEGER          , OPTIONAL ::  field_size19
REAL             , OPTIONAL ::  field19(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check20
INTEGER          , OPTIONAL ::  field_size20
REAL             , OPTIONAL ::  field20(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check21
INTEGER          , OPTIONAL ::  field_size21
REAL             , OPTIONAL ::  field21(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check22
INTEGER          , OPTIONAL ::  field_size22
REAL             , OPTIONAL ::  field22(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check23
INTEGER          , OPTIONAL ::  field_size23
REAL             , OPTIONAL ::  field23(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check24
INTEGER          , OPTIONAL ::  field_size24
REAL             , OPTIONAL ::  field24(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check25
INTEGER          , OPTIONAL ::  field_size25
REAL             , OPTIONAL ::  field25(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check26
INTEGER          , OPTIONAL ::  field_size26
REAL             , OPTIONAL ::  field26(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check27
INTEGER          , OPTIONAL ::  field_size27
REAL             , OPTIONAL ::  field27(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check28
INTEGER          , OPTIONAL ::  field_size28
REAL             , OPTIONAL ::  field28(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check29
INTEGER          , OPTIONAL ::  field_size29
REAL             , OPTIONAL ::  field29(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check30
INTEGER          , OPTIONAL ::  field_size30
REAL             , OPTIONAL ::  field30(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check31
INTEGER          , OPTIONAL ::  field_size31
REAL             , OPTIONAL ::  field31(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check32
INTEGER          , OPTIONAL ::  field_size32
REAL             , OPTIONAL ::  field32(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check33
INTEGER          , OPTIONAL ::  field_size33
REAL             , OPTIONAL ::  field33(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check34
INTEGER          , OPTIONAL ::  field_size34
REAL             , OPTIONAL ::  field34(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check35
INTEGER          , OPTIONAL ::  field_size35
REAL             , OPTIONAL ::  field35(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check36
INTEGER          , OPTIONAL ::  field_size36
REAL             , OPTIONAL ::  field36(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check37
INTEGER          , OPTIONAL ::  field_size37
REAL             , OPTIONAL ::  field37(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check38
INTEGER          , OPTIONAL ::  field_size38
REAL             , OPTIONAL ::  field38(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check39
INTEGER          , OPTIONAL ::  field_size39
REAL             , OPTIONAL ::  field39(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check40
INTEGER          , OPTIONAL ::  field_size40
REAL             , OPTIONAL ::  field40(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check41
INTEGER          , OPTIONAL ::  field_size41
REAL             , OPTIONAL ::  field41(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check42
INTEGER          , OPTIONAL ::  field_size42
REAL             , OPTIONAL ::  field42(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check43
INTEGER          , OPTIONAL ::  field_size43
REAL             , OPTIONAL ::  field43(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check44
INTEGER          , OPTIONAL ::  field_size44
REAL             , OPTIONAL ::  field44(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check45
INTEGER          , OPTIONAL ::  field_size45
REAL             , OPTIONAL ::  field45(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check46
INTEGER          , OPTIONAL ::  field_size46
REAL             , OPTIONAL ::  field46(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check47
INTEGER          , OPTIONAL ::  field_size47
REAL             , OPTIONAL ::  field47(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check48
INTEGER          , OPTIONAL ::  field_size48
REAL             , OPTIONAL ::  field48(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check49
INTEGER          , OPTIONAL ::  field_size49
REAL             , OPTIONAL ::  field49(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check50
INTEGER          , OPTIONAL ::  field_size50
REAL             , OPTIONAL ::  field50(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check51
INTEGER          , OPTIONAL ::  field_size51
REAL             , OPTIONAL ::  field51(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check52
INTEGER          , OPTIONAL ::  field_size52
REAL             , OPTIONAL ::  field52(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check53
INTEGER          , OPTIONAL ::  field_size53
REAL             , OPTIONAL ::  field53(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check54
INTEGER          , OPTIONAL ::  field_size54
REAL             , OPTIONAL ::  field54(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check55
INTEGER          , OPTIONAL ::  field_size55
REAL             , OPTIONAL ::  field55(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check56
INTEGER          , OPTIONAL ::  field_size56
REAL             , OPTIONAL ::  field56(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check57
INTEGER          , OPTIONAL ::  field_size57
REAL             , OPTIONAL ::  field57(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check58
INTEGER          , OPTIONAL ::  field_size58
REAL             , OPTIONAL ::  field58(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check59
INTEGER          , OPTIONAL ::  field_size59
REAL             , OPTIONAL ::  field59(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check60
INTEGER          , OPTIONAL ::  field_size60
REAL             , OPTIONAL ::  field60(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check61
INTEGER          , OPTIONAL ::  field_size61
REAL             , OPTIONAL ::  field61(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check62
INTEGER          , OPTIONAL ::  field_size62
REAL             , OPTIONAL ::  field62(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check63
INTEGER          , OPTIONAL ::  field_size63
REAL             , OPTIONAL ::  field63(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check64
INTEGER          , OPTIONAL ::  field_size64
REAL             , OPTIONAL ::  field64(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check65
INTEGER          , OPTIONAL ::  field_size65
REAL             , OPTIONAL ::  field65(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check66
INTEGER          , OPTIONAL ::  field_size66
REAL             , OPTIONAL ::  field66(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check67
INTEGER          , OPTIONAL ::  field_size67
REAL             , OPTIONAL ::  field67(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check68
INTEGER          , OPTIONAL ::  field_size68
REAL             , OPTIONAL ::  field68(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check69
INTEGER          , OPTIONAL ::  field_size69
REAL             , OPTIONAL ::  field69(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check70
INTEGER          , OPTIONAL ::  field_size70
REAL             , OPTIONAL ::  field70(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check71
INTEGER          , OPTIONAL ::  field_size71
REAL             , OPTIONAL ::  field71(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check72
INTEGER          , OPTIONAL ::  field_size72
REAL             , OPTIONAL ::  field72(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check73
INTEGER          , OPTIONAL ::  field_size73
REAL             , OPTIONAL ::  field73(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check74
INTEGER          , OPTIONAL ::  field_size74
REAL             , OPTIONAL ::  field74(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check75
INTEGER          , OPTIONAL ::  field_size75
REAL             , OPTIONAL ::  field75(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check76
INTEGER          , OPTIONAL ::  field_size76
REAL             , OPTIONAL ::  field76(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check77
INTEGER          , OPTIONAL ::  field_size77
REAL             , OPTIONAL ::  field77(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check78
INTEGER          , OPTIONAL ::  field_size78
REAL             , OPTIONAL ::  field78(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check79
INTEGER          , OPTIONAL ::  field_size79
REAL             , OPTIONAL ::  field79(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check80
INTEGER          , OPTIONAL ::  field_size80
REAL             , OPTIONAL ::  field80(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check81
INTEGER          , OPTIONAL ::  field_size81
REAL             , OPTIONAL ::  field81(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check82
INTEGER          , OPTIONAL ::  field_size82
REAL             , OPTIONAL ::  field82(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check83
INTEGER          , OPTIONAL ::  field_size83
REAL             , OPTIONAL ::  field83(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check84
INTEGER          , OPTIONAL ::  field_size84
REAL             , OPTIONAL ::  field84(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check85
INTEGER          , OPTIONAL ::  field_size85
REAL             , OPTIONAL ::  field85(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check86
INTEGER          , OPTIONAL ::  field_size86
REAL             , OPTIONAL ::  field86(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check87
INTEGER          , OPTIONAL ::  field_size87
REAL             , OPTIONAL ::  field87(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check88
INTEGER          , OPTIONAL ::  field_size88
REAL             , OPTIONAL ::  field88(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check89
INTEGER          , OPTIONAL ::  field_size89
REAL             , OPTIONAL ::  field89(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check90
INTEGER          , OPTIONAL ::  field_size90
REAL             , OPTIONAL ::  field90(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check91
INTEGER          , OPTIONAL ::  field_size91
REAL             , OPTIONAL ::  field91(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check92
INTEGER          , OPTIONAL ::  field_size92
REAL             , OPTIONAL ::  field92(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check93
INTEGER          , OPTIONAL ::  field_size93
REAL             , OPTIONAL ::  field93(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check94
INTEGER          , OPTIONAL ::  field_size94
REAL             , OPTIONAL ::  field94(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check95
INTEGER          , OPTIONAL ::  field_size95
REAL             , OPTIONAL ::  field95(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check96
INTEGER          , OPTIONAL ::  field_size96
REAL             , OPTIONAL ::  field96(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check97
INTEGER          , OPTIONAL ::  field_size97
REAL             , OPTIONAL ::  field97(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check98
INTEGER          , OPTIONAL ::  field_size98
REAL             , OPTIONAL ::  field98(*)
CHARACTER (len=5), OPTIONAL ::  key_to_check99
INTEGER          , OPTIONAL ::  field_size99
REAL             , OPTIONAL ::  field99(*)

IF (lhook) CALL dr_hook('INTEGRITY:CHECK_HASH_M',zhook_in,zhook_handle)

CALL check_hash(field ,field_size ,key_to_check)

IF(PRESENT(field1)) CALL check_hash(field1,field_size1,key_to_check1)
IF(PRESENT(field2)) CALL check_hash(field2,field_size2,key_to_check2)
IF(PRESENT(field3)) CALL check_hash(field3,field_size3,key_to_check3)
IF(PRESENT(field4)) CALL check_hash(field4,field_size4,key_to_check4)
IF(PRESENT(field5)) CALL check_hash(field5,field_size5,key_to_check5)
IF(PRESENT(field6)) CALL check_hash(field6,field_size6,key_to_check6)
IF(PRESENT(field7)) CALL check_hash(field7,field_size7,key_to_check7)
IF(PRESENT(field8)) CALL check_hash(field8,field_size8,key_to_check8)
IF(PRESENT(field9)) CALL check_hash(field9,field_size9,key_to_check9)
IF(PRESENT(field10)) CALL check_hash(field10,field_size10,key_to_check10)
IF(PRESENT(field11)) CALL check_hash(field11,field_size11,key_to_check11)
IF(PRESENT(field12)) CALL check_hash(field12,field_size12,key_to_check12)
IF(PRESENT(field13)) CALL check_hash(field13,field_size13,key_to_check13)
IF(PRESENT(field14)) CALL check_hash(field14,field_size14,key_to_check14)
IF(PRESENT(field15)) CALL check_hash(field15,field_size15,key_to_check15)
IF(PRESENT(field16)) CALL check_hash(field16,field_size16,key_to_check16)
IF(PRESENT(field17)) CALL check_hash(field17,field_size17,key_to_check17)
IF(PRESENT(field18)) CALL check_hash(field18,field_size18,key_to_check18)
IF(PRESENT(field19)) CALL check_hash(field19,field_size19,key_to_check19)
IF(PRESENT(field20)) CALL check_hash(field20,field_size20,key_to_check20)
IF(PRESENT(field21)) CALL check_hash(field21,field_size21,key_to_check21)
IF(PRESENT(field22)) CALL check_hash(field22,field_size22,key_to_check22)
IF(PRESENT(field23)) CALL check_hash(field23,field_size23,key_to_check23)
IF(PRESENT(field24)) CALL check_hash(field24,field_size24,key_to_check24)
IF(PRESENT(field25)) CALL check_hash(field25,field_size25,key_to_check25)
IF(PRESENT(field26)) CALL check_hash(field26,field_size26,key_to_check26)
IF(PRESENT(field27)) CALL check_hash(field27,field_size27,key_to_check27)
IF(PRESENT(field28)) CALL check_hash(field28,field_size28,key_to_check28)
IF(PRESENT(field29)) CALL check_hash(field29,field_size29,key_to_check29)
IF(PRESENT(field30)) CALL check_hash(field30,field_size30,key_to_check30)
IF(PRESENT(field31)) CALL check_hash(field31,field_size31,key_to_check31)
IF(PRESENT(field32)) CALL check_hash(field32,field_size32,key_to_check32)
IF(PRESENT(field33)) CALL check_hash(field33,field_size33,key_to_check33)
IF(PRESENT(field34)) CALL check_hash(field34,field_size34,key_to_check34)
IF(PRESENT(field35)) CALL check_hash(field35,field_size35,key_to_check35)
IF(PRESENT(field36)) CALL check_hash(field36,field_size36,key_to_check36)
IF(PRESENT(field37)) CALL check_hash(field37,field_size37,key_to_check37)
IF(PRESENT(field38)) CALL check_hash(field38,field_size38,key_to_check38)
IF(PRESENT(field39)) CALL check_hash(field39,field_size39,key_to_check39)
IF(PRESENT(field40)) CALL check_hash(field40,field_size40,key_to_check40)
IF(PRESENT(field41)) CALL check_hash(field41,field_size41,key_to_check41)
IF(PRESENT(field42)) CALL check_hash(field42,field_size42,key_to_check42)
IF(PRESENT(field43)) CALL check_hash(field43,field_size43,key_to_check43)
IF(PRESENT(field44)) CALL check_hash(field44,field_size44,key_to_check44)
IF(PRESENT(field45)) CALL check_hash(field45,field_size45,key_to_check45)
IF(PRESENT(field46)) CALL check_hash(field46,field_size46,key_to_check46)
IF(PRESENT(field47)) CALL check_hash(field47,field_size47,key_to_check47)
IF(PRESENT(field48)) CALL check_hash(field48,field_size48,key_to_check48)
IF(PRESENT(field49)) CALL check_hash(field49,field_size49,key_to_check49)
IF(PRESENT(field50)) CALL check_hash(field50,field_size50,key_to_check50)
IF(PRESENT(field51)) CALL check_hash(field51,field_size51,key_to_check51)
IF(PRESENT(field52)) CALL check_hash(field52,field_size52,key_to_check52)
IF(PRESENT(field53)) CALL check_hash(field53,field_size53,key_to_check53)
IF(PRESENT(field54)) CALL check_hash(field54,field_size54,key_to_check54)
IF(PRESENT(field55)) CALL check_hash(field55,field_size55,key_to_check55)
IF(PRESENT(field56)) CALL check_hash(field56,field_size56,key_to_check56)
IF(PRESENT(field57)) CALL check_hash(field57,field_size57,key_to_check57)
IF(PRESENT(field58)) CALL check_hash(field58,field_size58,key_to_check58)
IF(PRESENT(field59)) CALL check_hash(field59,field_size59,key_to_check59)
IF(PRESENT(field60)) CALL check_hash(field60,field_size60,key_to_check60)
IF(PRESENT(field61)) CALL check_hash(field61,field_size61,key_to_check61)
IF(PRESENT(field62)) CALL check_hash(field62,field_size62,key_to_check62)
IF(PRESENT(field63)) CALL check_hash(field63,field_size63,key_to_check63)
IF(PRESENT(field64)) CALL check_hash(field64,field_size64,key_to_check64)
IF(PRESENT(field65)) CALL check_hash(field65,field_size65,key_to_check65)
IF(PRESENT(field66)) CALL check_hash(field66,field_size66,key_to_check66)
IF(PRESENT(field67)) CALL check_hash(field67,field_size67,key_to_check67)
IF(PRESENT(field68)) CALL check_hash(field68,field_size68,key_to_check68)
IF(PRESENT(field69)) CALL check_hash(field69,field_size69,key_to_check69)
IF(PRESENT(field70)) CALL check_hash(field70,field_size70,key_to_check70)
IF(PRESENT(field71)) CALL check_hash(field71,field_size71,key_to_check71)
IF(PRESENT(field72)) CALL check_hash(field72,field_size72,key_to_check72)
IF(PRESENT(field73)) CALL check_hash(field73,field_size73,key_to_check73)
IF(PRESENT(field74)) CALL check_hash(field74,field_size74,key_to_check74)
IF(PRESENT(field75)) CALL check_hash(field75,field_size75,key_to_check75)
IF(PRESENT(field76)) CALL check_hash(field76,field_size76,key_to_check76)
IF(PRESENT(field77)) CALL check_hash(field77,field_size77,key_to_check77)
IF(PRESENT(field78)) CALL check_hash(field78,field_size78,key_to_check78)
IF(PRESENT(field79)) CALL check_hash(field79,field_size79,key_to_check79)
IF(PRESENT(field80)) CALL check_hash(field80,field_size80,key_to_check80)
IF(PRESENT(field81)) CALL check_hash(field81,field_size81,key_to_check81)
IF(PRESENT(field82)) CALL check_hash(field82,field_size82,key_to_check82)
IF(PRESENT(field83)) CALL check_hash(field83,field_size83,key_to_check83)
IF(PRESENT(field84)) CALL check_hash(field84,field_size84,key_to_check84)
IF(PRESENT(field85)) CALL check_hash(field85,field_size85,key_to_check85)
IF(PRESENT(field86)) CALL check_hash(field86,field_size86,key_to_check86)
IF(PRESENT(field87)) CALL check_hash(field87,field_size87,key_to_check87)
IF(PRESENT(field88)) CALL check_hash(field88,field_size88,key_to_check88)
IF(PRESENT(field89)) CALL check_hash(field89,field_size89,key_to_check89)
IF(PRESENT(field90)) CALL check_hash(field90,field_size90,key_to_check90)
IF(PRESENT(field91)) CALL check_hash(field91,field_size91,key_to_check91)
IF(PRESENT(field92)) CALL check_hash(field92,field_size92,key_to_check92)
IF(PRESENT(field93)) CALL check_hash(field93,field_size93,key_to_check93)
IF(PRESENT(field94)) CALL check_hash(field94,field_size94,key_to_check94)
IF(PRESENT(field95)) CALL check_hash(field95,field_size95,key_to_check95)
IF(PRESENT(field96)) CALL check_hash(field96,field_size96,key_to_check96)
IF(PRESENT(field97)) CALL check_hash(field97,field_size97,key_to_check97)
IF(PRESENT(field98)) CALL check_hash(field98,field_size98,key_to_check98)
IF(PRESENT(field99)) CALL check_hash(field99,field_size99,key_to_check99)

IF (lhook) CALL dr_hook('INTEGRITY:CHECK_HASH_M',zhook_out,zhook_handle)

END SUBROUTINE

SUBROUTINE add_hash_m(field ,field_size ,key_to_add,&
                      field1,field_size1,key_to_add1,&
                      field2,field_size2,key_to_add2,&
                      field3,field_size3,key_to_add3,&
                      field4,field_size4,key_to_add4,&
                      field5,field_size5,key_to_add5,&
                      field6,field_size6,key_to_add6,&
                      field7,field_size7,key_to_add7,&
                      field8,field_size8,key_to_add8,&
                      field9,field_size9,key_to_add9)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (len=5)  key_to_add
INTEGER            field_size
REAL               field(*)

CHARACTER (len=5), OPTIONAL ::  key_to_add1
INTEGER          , OPTIONAL ::  field_size1
REAL             , OPTIONAL ::  field1(*)
CHARACTER (len=5), OPTIONAL ::  key_to_add2
INTEGER          , OPTIONAL ::  field_size2
REAL             , OPTIONAL ::  field2(*)
CHARACTER (len=5), OPTIONAL ::  key_to_add3
INTEGER          , OPTIONAL ::  field_size3
REAL             , OPTIONAL ::  field3(*)
CHARACTER (len=5), OPTIONAL ::  key_to_add4
INTEGER          , OPTIONAL ::  field_size4
REAL             , OPTIONAL ::  field4(*)
CHARACTER (len=5), OPTIONAL ::  key_to_add5
INTEGER          , OPTIONAL ::  field_size5
REAL             , OPTIONAL ::  field5(*)
CHARACTER (len=5), OPTIONAL ::  key_to_add6
INTEGER          , OPTIONAL ::  field_size6
REAL             , OPTIONAL ::  field6(*)
CHARACTER (len=5), OPTIONAL ::  key_to_add7
INTEGER          , OPTIONAL ::  field_size7
REAL             , OPTIONAL ::  field7(*)
CHARACTER (len=5), OPTIONAL ::  key_to_add8
INTEGER          , OPTIONAL ::  field_size8
REAL             , OPTIONAL ::  field8(*)
CHARACTER (len=5), OPTIONAL ::  key_to_add9
INTEGER          , OPTIONAL ::  field_size9
REAL             , OPTIONAL ::  field9(*)

IF (lhook) CALL dr_hook('INTEGRITY:ADD_HASH_M',zhook_in,zhook_handle)

CALL add_hash(field ,field_size ,key_to_add)

IF(PRESENT(field1)) CALL add_hash(field1,field_size1,key_to_add1)
IF(PRESENT(field2)) CALL add_hash(field2,field_size2,key_to_add2)
IF(PRESENT(field3)) CALL add_hash(field3,field_size3,key_to_add3)
IF(PRESENT(field4)) CALL add_hash(field4,field_size4,key_to_add4)
IF(PRESENT(field5)) CALL add_hash(field5,field_size5,key_to_add5)
IF(PRESENT(field6)) CALL add_hash(field6,field_size6,key_to_add6)
IF(PRESENT(field7)) CALL add_hash(field7,field_size7,key_to_add7)
IF(PRESENT(field8)) CALL add_hash(field8,field_size8,key_to_add8)
IF(PRESENT(field9)) CALL add_hash(field9,field_size9,key_to_add9)

IF (lhook) CALL dr_hook('INTEGRITY:ADD_HASH_M',zhook_out,zhook_handle)

END SUBROUTINE


SUBROUTINE reset_AV()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE proc_info_mod, ONLY : me

IMPLICIT NONE 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER i

IF (lhook) CALL dr_hook('INTEGRITY:RESET_AV',zhook_in,zhook_handle)

last_entry=0
! initialise (in case they key is less than 65 characters,as
! on our old 1 system)
 DO i=1,registry_length+1
   hash(i) = "                                                                  "
 END DO

 DO i=1,registry_length
   key(i) = "     "
 END DO

IF (me == 0) WRITE(av_out,fmt='(A)') 'RESETTING AV ++++++++++++++++++'
IF (lhook) CALL dr_hook('INTEGRITY:RESET_AV',zhook_out,zhook_handle)

END SUBROUTINE reset_AV


SUBROUTINE add_hash(field,field_size,key_to_add)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE proc_info_mod,     ONLY : me,n_proc
USE ereport_mod, ONLY : ereport
USE timestep_mod,  ONLY : timestep_number

IMPLICIT NONE 
CHARACTER (len=9) filename
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (len=5)  key_to_add
INTEGER field_size,i,ierr
REAL    field(*),field_max,field_min

CHARACTER(len=80) routinename, cmessage
CHARACTER(len=65) hashvalue 
!
!  Purpose: Defines unit numbers relevant to history file
!           and variables used to hold the logical to physical
!           file associations made within the model
!
!  Logical Filenames used in the model
!
      CHARACTER(LEN=256) hkfile,ppxref,config,stashctl,namelist,output,      &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,ftxx,    &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar
     

!
      CHARACTER(LEN=256) MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                                 ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        RSUB_UNIT,                                                &
                                 ! File indicating whether resub required
     &        XHIST_UNIT,                                               &
                                 ! Main history file unit
     &        THIST_UNIT,                                               &
                                 ! Backup history file unit
     &        HKFILE_UNIT,                                              &
                                 ! Operational houskeeping file unit    
     &        EG_UNIT            ! ENDGame diagnostics/info unit
!
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(RSUB_UNIT =10)
      PARAMETER(XHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)

!
! Parameters specifying unit numbers relevant to ENDGame diagnostics
!
      PARAMETER(EG_UNIT  = 55)

! UKCA unit numbers

      INTEGER, PARAMETER :: ukcafjxx_unit=170 ! Fast-J(X) cross section data
      INTEGER, PARAMETER :: ukcafjsc_unit=171 ! Fast-JX scattering data
      INTEGER, PARAMETER :: ukca2do3_unit=172 ! 2D top boundary O3 data 
      INTEGER, PARAMETER :: ukca2ch4_unit=173 ! 2D top boundary CH4 data
      INTEGER, PARAMETER :: ukca2noy_unit=174 ! 2D top boundary NOY data
      INTEGER, PARAMETER :: ukca2pho_unit=175 ! 2D photolysis input data
      INTEGER, PARAMETER :: ukcastrd_unit=176 ! Stratospheric model radiation field. 
      INTEGER, PARAMETER :: ukcasto3_unit=177 ! Strat standard atmosphere T and O3.
      INTEGER, PARAMETER :: ukcastar_unit=178 ! Stratospheric sulfate aerosol climatology 
      INTEGER, PARAMETER :: ukcafjar_unit=179 ! Sulfate aerosol cliamtology for Fast-JX
! Text output file for STASH-related information is assigned to UNIT 200

!
! Namelist of all permissible logical files.
!
      NAMELIST / nlcfiles /                                             &
                   hkfile,ppxref,config,stashctl,namelist,output,       &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,         &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(RSUB       ,MODEL_FT_UNIT(10) ), &
     &(XHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(ICECALVE  ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &                                (ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(GENLAND   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
                                      (lwspectd   ,model_ft_unit(80) ), &
      (surgeou1  ,model_ft_unit(81) ),(surgeout   ,model_ft_unit(82) ), &
      (ppscreen  ,model_ft_unit(83) ),(ppsmc      ,model_ft_unit(84) ), &
      (wfout     ,model_ft_unit(85) ),(uarsout1   ,model_ft_unit(86) ), &
      (uarsout2  ,model_ft_unit(87) ),(icefout    ,model_ft_unit(88) ), &
      (mosout    ,model_ft_unit(89) ),(vert_lev   ,model_ft_unit(90) ), &
      (sstout    ,model_ft_unit(91) ),(siceout    ,model_ft_unit(92) ), &
      (curntout  ,model_ft_unit(93) ),(flxcrout   ,model_ft_unit(94) ), &
      (dmsconc   ,model_ft_unit(95) ),(orog       ,model_ft_unit(96) ), &
      (transp    ,model_ft_unit(97) ),(olabcin    ,model_ft_unit(98) ), &
      (ocndepth  ,model_ft_unit(99) ),                                  &
      (foamout1  ,model_ft_unit(100)),(foamout2   ,model_ft_unit(101)), &
      (cxbkgerr  ,model_ft_unit(102)),(rfmout     ,model_ft_unit(103)), &
      (idealise  ,model_ft_unit(106)),(tdf_dump   ,model_ft_unit(107)), &
      (iau_inc   ,model_ft_unit(108)),(murkfile   ,model_ft_unit(109)), &
      (sulpemis  ,model_ft_unit(110)),(usrancil   ,model_ft_unit(111)), &
      (usrmulti  ,model_ft_unit(112)),(ousrancl   ,model_ft_unit(113)), &
      (ousrmult  ,model_ft_unit(114)),(so2natem   ,model_ft_unit(115)), &
      (chemoxid  ,model_ft_unit(116)),(aerofcg    ,model_ft_unit(117)), &
      (co2emits  ,model_ft_unit(118)),(tppsozon   ,model_ft_unit(119)), &
      (landfrac  ,model_ft_unit(120)),(wlabcou1   ,model_ft_unit(121)), &
      (wlabcou2  ,model_ft_unit(122)),(wlabcou3   ,model_ft_unit(123)), &
      (wlabcou4  ,model_ft_unit(124)),(alabcin1   ,model_ft_unit(125)), &
      (alabcin2  ,model_ft_unit(126)),                                  &
      (ocffemis  ,model_ft_unit(128)),(horzgrid   ,model_ft_unit(129)), &
      (surfemis  ,model_ft_unit(130)),(aircrems   ,model_ft_unit(131)), &
      (stratems  ,model_ft_unit(132)),(extraems   ,model_ft_unit(133)), &
      (radonems  ,model_ft_unit(134)),(fracinit   ,model_ft_unit(135)), &
      (veginit   ,model_ft_unit(136)),(disturb    ,model_ft_unit(137)), &
      (cached    ,model_ft_unit(138)),(sootemis   ,model_ft_unit(139)), &
      (alabcou1  ,model_ft_unit(140)),(alabcou2   ,model_ft_unit(141)), &
      (alabcou3  ,model_ft_unit(142)),(alabcou4   ,model_ft_unit(143)), &
      (alabcou5  ,model_ft_unit(144)),(alabcou6   ,model_ft_unit(145)), &
      (alabcou7  ,model_ft_unit(146)),(alabcou8   ,model_ft_unit(147)), &
      (cariolo3  ,model_ft_unit(148)),(rpseed     ,model_ft_unit(149)), &
      (ppvar     ,model_ft_unit(150)),(pp10       ,model_ft_unit(151)), &
      (icfile    ,model_ft_unit(152)),(var_grid   ,model_ft_unit(153)), &
      (arclbiog  ,model_ft_unit(154)),(arclbiom   ,model_ft_unit(155)), &
      (arclblck  ,model_ft_unit(156)),(arclsslt   ,model_ft_unit(157)), &
      (arclsulp  ,model_ft_unit(158)),(arcldust   ,model_ft_unit(159)), &
      (arclocff  ,model_ft_unit(160)),(arcldlta   ,model_ft_unit(161)), &
      (topmean   ,model_ft_unit(162)),(topstdev   ,model_ft_unit(163)), &
      (ppmbc     ,model_ft_unit(164)),(ukcaprec   ,model_ft_unit(165)), &
      (ukcaacsw  ,model_ft_unit(166)),(ukcaaclw   ,model_ft_unit(167)), &
      (ukcacrsw  ,model_ft_unit(168)),(ukcacrlw   ,model_ft_unit(169)), &
      (ukcafjxx  ,model_ft_unit(170)),(ukcafjsc   ,model_ft_unit(171)), &
      (ukca2do3  ,model_ft_unit(172)),(ukca2ch4   ,model_ft_unit(173)), &
      (ukca2noy  ,model_ft_unit(174)),(ukca2pho   ,model_ft_unit(175)), &
      (ukcastrd  ,model_ft_unit(176)),(ukcasto3   ,model_ft_unit(177)), &
      (ukcastar  ,model_ft_unit(178)),(ukcafjar   ,model_ft_unit(179))
! Text output file for STASH-related information is assigned to UNIT 200


IF (lhook) CALL dr_hook('INTEGRITY:ADD_HASH',zhook_in,zhook_handle)

routinename='add_hash (integrity_mod)'

! check for a NaN in the field to add

IF(ANY(field(1:field_size)/=field(1:field_size))) THEN
 

  IF (me == 0) WRITE(av_out,fmt='(2A,I20)') key_to_add,                    &
                       ': NaN detected! Field size on add :',field_size

!   cmessage = ' NaN detected in add of ' // TRIM(key_to_add)
!   CALL ereport(routinename,1,cmessage)

ENDIF

IF(last_entry==0) THEN
! initialise (in case they key is less than 65 characters)
 DO i=1,registry_length+1
   hash(i) = "                                                                  "
 END DO

ENDIF

DO i=1,last_entry
  IF(key(i) == key_to_add) THEN
   WRITE(av_out,*) 'ATTEMPTING TO ADD DUPLICATE KEY!',key_to_add
   cmessage = 'Attempted to add duplicate key'
   CALL ereport(RoutineName, 1, cmessage)
  END IF
END DO

last_entry = last_entry + 1

IF (last_entry.le.registry_length) THEN

  key(last_entry) = key_to_add
  ! DEPENDS ON: eg_hash
  CALL eg_hash(field,8*field_size,hash(last_entry))

  field_max=maxval(field(1:field_size))
  field_min=minval(field(1:field_size))

  CALL gc_rmax(1,n_proc,ierr,field_max)
  CALL gc_rmin(1,n_proc,ierr,field_min)

  hashvalue = hash(last_entry) 

  IF (me == 0) WRITE(av_out,fmt='(A,I4.4,A,A,I5,A,I5.5,2E32.16,2A)')  &
                            'adding on PE '                           &
                            ,me,' : ', key_to_add,last_entry,         &
                          ' ',field_size,field_min,                   &
               field_max,' ',hashvalue(1:64)
  IF (integrity_file) THEN
    WRITE (filename,fmt='(I4.4,A,I4.4)') me,'-',timestep_number
    OPEN(eg_unit,file=filename, position='APPEND')
    WRITE(eg_unit,'(A,A)') key_to_add,hashvalue(1:63)
    CLOSE(eg_unit)
  END IF
  !  
  ! regeneratate hash table verify hash
  ! DEPENDS ON: eg_hash
  CALL eg_hash(hash,65*(last_entry),hash(registry_length+1))

ELSE
  WRITE(0,*) 'REGISTRY FULL!'
END IF

IF (lhook) CALL dr_hook('INTEGRITY:ADD_HASH',zhook_out,zhook_handle)

END SUBROUTINE


SUBROUTINE check_hash(field,field_size,key_to_check,debug_text)

! this #if is needed because the module needs to be
! excluded in the configuration (when in use)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE proc_info_mod,     ONLY : me
USE ereport_mod, ONLY : ereport

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (len=5)  key_to_check
INTEGER field_size,entry
REAL    field(*)
CHARACTER(len=65) hash_to_check
CHARACTER(LEN=*) , OPTIONAL :: debug_text

CHARACTER(len=80) routinename,cmessage


! this #if is needed because the include file needs to be 
! excluded in the configuration (when in use)

IF (lhook) CALL dr_hook('INTEGRITY:CHECK_HASH',zhook_in,zhook_handle)

routinename='check_hash (integrity_mod)'

! check for a NaN in the field to check:

IF(ANY(field(1:field_size).ne.field(1:field_size))) THEN

   IF (me == 0) WRITE(av_out,fmt='(2A,I20)') key_to_check,&
                  ': NaN detected! Field size on check :',field_size

!
! do not abort on NaN detection at the moment,
! as NaN's are picked up at the beginning where
! fields have not yet been propperly initialised
!
! #if required as these calls are otherwise not resolved

!#if defined (IBM) && (EG_INTEGRITY_TEST)
!        call xl__trbk
!#endif
!#if defined (1) && (EG_INTEGRITY_TEST)
!CALL TRACEBACKQQ(user_exit_code=0)
!#endif

!   cmessage = ' NaN detected in check of ' // TRIM(key_to_check)

!   CALL ereport(routinename,1,cmessage)

ENDIF

hash_to_check=&
   "                                                                  "

entry=1
DO while (key(entry) /= key_to_check .AND. entry.lt.registry_length)
  entry=entry+1
END DO

IF(entry.le.registry_length.and.key(entry) == key_to_check) THEN

  ! check registry itself first
  ! DEPENDS ON: eg_hash
  CALL eg_hash(hash,65*(last_entry),hash_to_check)
  IF(hash_to_check.ne.hash(registry_length+1)) THEN
   IF (me == 0) WRITE(0,fmt='(A)') 'registry corrupted'
!
! dump a full trace back and abort (so that we know where things went wrong!)
! #if required as we do not link against the correct libraries by default


  END IF

  ! DEPENDS ON: eg_hash
  CALL eg_hash(field,8*field_size,hash_to_check)

  IF(PRESENT(debug_text))  WRITE(av_out,fmt='(2A)',ADVANCE='NO') debug_text,' '

  IF(hash_to_check==hash(entry)) THEN
    IF (me == 0) THEN

      WRITE(av_out,fmt='(A,I4.4,3A,I20,2A)') 'PE ',me,' ', key_to_check, &
              ': ok! Field size/ hash on check :',                       &
              field_size,' ',hash_to_check(1:64)
    END IF
  ELSE
   WRITE(av_out,fmt='(2A,I20)') key_to_check,                            &
         ': corrupted! Field size on check :',field_size

! #if required as we do not link against the correct libraries by default


   cmessage = '      '
 
   IF(PRESENT(debug_text)) cmessage = TRIM(debug_text) // ' '

   cmessage = TRIM(cmessage) // TRIM(key_to_check) // ': corrupted'

   CALL ereport(routinename,field_size,cmessage)

  END IF
ELSE
  IF (me == 0) WRITE(av_out,fmt='(3A)') 'key ',key_to_check,' not found'
END IF

IF (lhook) CALL dr_hook('INTEGRITY:CHECK_HASH',zhook_out,zhook_handle)

END SUBROUTINE

SUBROUTINE update_hash_m(field,field_size,key_to_update,&
                       field1,field_size1,key_to_update1,&
                       field2,field_size2,key_to_update2,&
                       field3,field_size3,key_to_update3,&
                       field4,field_size4,key_to_update4,&
                       field5,field_size5,key_to_update5,&
                       field6,field_size6,key_to_update6,&
                       field7,field_size7,key_to_update7,&
                       field8,field_size8,key_to_update8,&
                       field9,field_size9,key_to_update9,&
                       field10 ,field_size10 ,key_to_update10,&
                       field11 ,field_size11 ,key_to_update11,&
                       field12 ,field_size12 ,key_to_update12,&
                       field13 ,field_size13 ,key_to_update13,&
                       field14 ,field_size14 ,key_to_update14,&
                       field15 ,field_size15 ,key_to_update15,&
                       field16 ,field_size16 ,key_to_update16,&
                       field17 ,field_size17 ,key_to_update17,&
                       field18 ,field_size18 ,key_to_update18,&
                       field19 ,field_size19 ,key_to_update19,&
                       field20 ,field_size20 ,key_to_update20,&
                       field21 ,field_size21 ,key_to_update21,&
                       field22 ,field_size22 ,key_to_update22,&
                       field23 ,field_size23 ,key_to_update23,&
                       field24 ,field_size24 ,key_to_update24,&
                       field25 ,field_size25 ,key_to_update25,&
                       field26 ,field_size26 ,key_to_update26,&
                       field27 ,field_size27 ,key_to_update27,&
                       field28 ,field_size28 ,key_to_update28,&
                       field29 ,field_size29 ,key_to_update29,&
                       field30 ,field_size30 ,key_to_update30,&
                       field31 ,field_size31 ,key_to_update31,&
                       field32 ,field_size32 ,key_to_update32,&
                       field33 ,field_size33 ,key_to_update33,&
                       field34 ,field_size34 ,key_to_update34,&
                       field35 ,field_size35 ,key_to_update35,&
                       field36 ,field_size36 ,key_to_update36,&
                       field37 ,field_size37 ,key_to_update37,&
                       field38 ,field_size38 ,key_to_update38,&
                       field39 ,field_size39 ,key_to_update39,&
                       field40 ,field_size40 ,key_to_update40,&
                       field41 ,field_size41 ,key_to_update41,&
                       field42 ,field_size42 ,key_to_update42,&
                       field43 ,field_size43 ,key_to_update43,&
                       field44 ,field_size44 ,key_to_update44,&
                       field45 ,field_size45 ,key_to_update45,&
                       field46 ,field_size46 ,key_to_update46,&
                       field47 ,field_size47 ,key_to_update47,&
                       field48 ,field_size48 ,key_to_update48,&
                       field49 ,field_size49 ,key_to_update49,&
                       field50 ,field_size50 ,key_to_update50,&
                       field51 ,field_size51 ,key_to_update51,&
                       field52 ,field_size52 ,key_to_update52,&
                       field53 ,field_size53 ,key_to_update53,&
                       field54 ,field_size54 ,key_to_update54,&
                       field55 ,field_size55 ,key_to_update55,&
                       field56 ,field_size56 ,key_to_update56,&
                       field57 ,field_size57 ,key_to_update57,&
                       field58 ,field_size58 ,key_to_update58,&
                       field59 ,field_size59 ,key_to_update59,&
                       field60 ,field_size60 ,key_to_update60,&
                       field61 ,field_size61 ,key_to_update61,&
                       field62 ,field_size62 ,key_to_update62,&
                       field63 ,field_size63 ,key_to_update63,&
                       field64 ,field_size64 ,key_to_update64,&
                       field65 ,field_size65 ,key_to_update65,&
                       field66 ,field_size66 ,key_to_update66,&
                       field67 ,field_size67 ,key_to_update67,&
                       field68 ,field_size68 ,key_to_update68,&
                       field69 ,field_size69 ,key_to_update69,&
                       field70 ,field_size70 ,key_to_update70,&
                       field71 ,field_size71 ,key_to_update71,&
                       field72 ,field_size72 ,key_to_update72,&
                       field73 ,field_size73 ,key_to_update73,&
                       field74 ,field_size74 ,key_to_update74,&
                       field75 ,field_size75 ,key_to_update75,&
                       field76 ,field_size76 ,key_to_update76,&
                       field77 ,field_size77 ,key_to_update77,&
                       field78 ,field_size78 ,key_to_update78,&
                       field79 ,field_size79 ,key_to_update79,&
                       field80 ,field_size80 ,key_to_update80,&
                       field81 ,field_size81 ,key_to_update81,&
                       field82 ,field_size82 ,key_to_update82,&
                       field83 ,field_size83 ,key_to_update83,&
                       field84 ,field_size84 ,key_to_update84,&
                       field85 ,field_size85 ,key_to_update85,&
                       field86 ,field_size86 ,key_to_update86,&
                       field87 ,field_size87 ,key_to_update87,&
                       field88 ,field_size88 ,key_to_update88,&
                       field89 ,field_size89 ,key_to_update89,&
                       field90 ,field_size90 ,key_to_update90,&
                       field91 ,field_size91 ,key_to_update91,&
                       field92 ,field_size92 ,key_to_update92,&
                       field93 ,field_size93 ,key_to_update93,&
                       field94 ,field_size94 ,key_to_update94,&
                       field95 ,field_size95 ,key_to_update95,&
                       field96 ,field_size96 ,key_to_update96,&
                       field97 ,field_size97 ,key_to_update97,&
                       field98 ,field_size98 ,key_to_update98,&
                       field99 ,field_size99 ,key_to_update99) 

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (len=5)  key_to_update
INTEGER            field_size
REAL               field(*)

CHARACTER (len=5), OPTIONAL ::  key_to_update1
INTEGER          , OPTIONAL ::  field_size1
REAL             , OPTIONAL ::  field1(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update2
INTEGER          , OPTIONAL ::  field_size2
REAL             , OPTIONAL ::  field2(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update3
INTEGER          , OPTIONAL ::  field_size3
REAL             , OPTIONAL ::  field3(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update4
INTEGER          , OPTIONAL ::  field_size4
REAL             , OPTIONAL ::  field4(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update5
INTEGER          , OPTIONAL ::  field_size5
REAL             , OPTIONAL ::  field5(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update6
INTEGER          , OPTIONAL ::  field_size6
REAL             , OPTIONAL ::  field6(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update7
INTEGER          , OPTIONAL ::  field_size7
REAL             , OPTIONAL ::  field7(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update8
INTEGER          , OPTIONAL ::  field_size8
REAL             , OPTIONAL ::  field8(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update9
INTEGER          , OPTIONAL ::  field_size9
REAL             , OPTIONAL ::  field9(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update10
INTEGER          , OPTIONAL ::  field_size10
REAL             , OPTIONAL ::  field10(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update11
INTEGER          , OPTIONAL ::  field_size11
REAL             , OPTIONAL ::  field11(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update12
INTEGER          , OPTIONAL ::  field_size12
REAL             , OPTIONAL ::  field12(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update13
INTEGER          , OPTIONAL ::  field_size13
REAL             , OPTIONAL ::  field13(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update14
INTEGER          , OPTIONAL ::  field_size14
REAL             , OPTIONAL ::  field14(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update15
INTEGER          , OPTIONAL ::  field_size15
REAL             , OPTIONAL ::  field15(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update16
INTEGER          , OPTIONAL ::  field_size16
REAL             , OPTIONAL ::  field16(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update17
INTEGER          , OPTIONAL ::  field_size17
REAL             , OPTIONAL ::  field17(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update18
INTEGER          , OPTIONAL ::  field_size18
REAL             , OPTIONAL ::  field18(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update19
INTEGER          , OPTIONAL ::  field_size19
REAL             , OPTIONAL ::  field19(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update20
INTEGER          , OPTIONAL ::  field_size20
REAL             , OPTIONAL ::  field20(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update21
INTEGER          , OPTIONAL ::  field_size21
REAL             , OPTIONAL ::  field21(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update22
INTEGER          , OPTIONAL ::  field_size22
REAL             , OPTIONAL ::  field22(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update23
INTEGER          , OPTIONAL ::  field_size23
REAL             , OPTIONAL ::  field23(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update24
INTEGER          , OPTIONAL ::  field_size24
REAL             , OPTIONAL ::  field24(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update25
INTEGER          , OPTIONAL ::  field_size25
REAL             , OPTIONAL ::  field25(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update26
INTEGER          , OPTIONAL ::  field_size26
REAL             , OPTIONAL ::  field26(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update27
INTEGER          , OPTIONAL ::  field_size27
REAL             , OPTIONAL ::  field27(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update28
INTEGER          , OPTIONAL ::  field_size28
REAL             , OPTIONAL ::  field28(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update29
INTEGER          , OPTIONAL ::  field_size29
REAL             , OPTIONAL ::  field29(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update30
INTEGER          , OPTIONAL ::  field_size30
REAL             , OPTIONAL ::  field30(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update31
INTEGER          , OPTIONAL ::  field_size31
REAL             , OPTIONAL ::  field31(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update32
INTEGER          , OPTIONAL ::  field_size32
REAL             , OPTIONAL ::  field32(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update33
INTEGER          , OPTIONAL ::  field_size33
REAL             , OPTIONAL ::  field33(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update34
INTEGER          , OPTIONAL ::  field_size34
REAL             , OPTIONAL ::  field34(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update35
INTEGER          , OPTIONAL ::  field_size35
REAL             , OPTIONAL ::  field35(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update36
INTEGER          , OPTIONAL ::  field_size36
REAL             , OPTIONAL ::  field36(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update37
INTEGER          , OPTIONAL ::  field_size37
REAL             , OPTIONAL ::  field37(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update38
INTEGER          , OPTIONAL ::  field_size38
REAL             , OPTIONAL ::  field38(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update39
INTEGER          , OPTIONAL ::  field_size39
REAL             , OPTIONAL ::  field39(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update40
INTEGER          , OPTIONAL ::  field_size40
REAL             , OPTIONAL ::  field40(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update41
INTEGER          , OPTIONAL ::  field_size41
REAL             , OPTIONAL ::  field41(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update42
INTEGER          , OPTIONAL ::  field_size42
REAL             , OPTIONAL ::  field42(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update43
INTEGER          , OPTIONAL ::  field_size43
REAL             , OPTIONAL ::  field43(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update44
INTEGER          , OPTIONAL ::  field_size44
REAL             , OPTIONAL ::  field44(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update45
INTEGER          , OPTIONAL ::  field_size45
REAL             , OPTIONAL ::  field45(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update46
INTEGER          , OPTIONAL ::  field_size46
REAL             , OPTIONAL ::  field46(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update47
INTEGER          , OPTIONAL ::  field_size47
REAL             , OPTIONAL ::  field47(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update48
INTEGER          , OPTIONAL ::  field_size48
REAL             , OPTIONAL ::  field48(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update49
INTEGER          , OPTIONAL ::  field_size49
REAL             , OPTIONAL ::  field49(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update50
INTEGER          , OPTIONAL ::  field_size50
REAL             , OPTIONAL ::  field50(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update51
INTEGER          , OPTIONAL ::  field_size51
REAL             , OPTIONAL ::  field51(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update52
INTEGER          , OPTIONAL ::  field_size52
REAL             , OPTIONAL ::  field52(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update53
INTEGER          , OPTIONAL ::  field_size53
REAL             , OPTIONAL ::  field53(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update54
INTEGER          , OPTIONAL ::  field_size54
REAL             , OPTIONAL ::  field54(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update55
INTEGER          , OPTIONAL ::  field_size55
REAL             , OPTIONAL ::  field55(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update56
INTEGER          , OPTIONAL ::  field_size56
REAL             , OPTIONAL ::  field56(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update57
INTEGER          , OPTIONAL ::  field_size57
REAL             , OPTIONAL ::  field57(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update58
INTEGER          , OPTIONAL ::  field_size58
REAL             , OPTIONAL ::  field58(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update59
INTEGER          , OPTIONAL ::  field_size59
REAL             , OPTIONAL ::  field59(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update60
INTEGER          , OPTIONAL ::  field_size60
REAL             , OPTIONAL ::  field60(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update61
INTEGER          , OPTIONAL ::  field_size61
REAL             , OPTIONAL ::  field61(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update62
INTEGER          , OPTIONAL ::  field_size62
REAL             , OPTIONAL ::  field62(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update63
INTEGER          , OPTIONAL ::  field_size63
REAL             , OPTIONAL ::  field63(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update64
INTEGER          , OPTIONAL ::  field_size64
REAL             , OPTIONAL ::  field64(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update65
INTEGER          , OPTIONAL ::  field_size65
REAL             , OPTIONAL ::  field65(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update66
INTEGER          , OPTIONAL ::  field_size66
REAL             , OPTIONAL ::  field66(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update67
INTEGER          , OPTIONAL ::  field_size67
REAL             , OPTIONAL ::  field67(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update68
INTEGER          , OPTIONAL ::  field_size68
REAL             , OPTIONAL ::  field68(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update69
INTEGER          , OPTIONAL ::  field_size69
REAL             , OPTIONAL ::  field69(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update70
INTEGER          , OPTIONAL ::  field_size70
REAL             , OPTIONAL ::  field70(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update71
INTEGER          , OPTIONAL ::  field_size71
REAL             , OPTIONAL ::  field71(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update72
INTEGER          , OPTIONAL ::  field_size72
REAL             , OPTIONAL ::  field72(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update73
INTEGER          , OPTIONAL ::  field_size73
REAL             , OPTIONAL ::  field73(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update74
INTEGER          , OPTIONAL ::  field_size74
REAL             , OPTIONAL ::  field74(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update75
INTEGER          , OPTIONAL ::  field_size75
REAL             , OPTIONAL ::  field75(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update76
INTEGER          , OPTIONAL ::  field_size76
REAL             , OPTIONAL ::  field76(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update77
INTEGER          , OPTIONAL ::  field_size77
REAL             , OPTIONAL ::  field77(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update78
INTEGER          , OPTIONAL ::  field_size78
REAL             , OPTIONAL ::  field78(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update79
INTEGER          , OPTIONAL ::  field_size79
REAL             , OPTIONAL ::  field79(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update80
INTEGER          , OPTIONAL ::  field_size80
REAL             , OPTIONAL ::  field80(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update81
INTEGER          , OPTIONAL ::  field_size81
REAL             , OPTIONAL ::  field81(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update82
INTEGER          , OPTIONAL ::  field_size82
REAL             , OPTIONAL ::  field82(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update83
INTEGER          , OPTIONAL ::  field_size83
REAL             , OPTIONAL ::  field83(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update84
INTEGER          , OPTIONAL ::  field_size84
REAL             , OPTIONAL ::  field84(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update85
INTEGER          , OPTIONAL ::  field_size85
REAL             , OPTIONAL ::  field85(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update86
INTEGER          , OPTIONAL ::  field_size86
REAL             , OPTIONAL ::  field86(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update87
INTEGER          , OPTIONAL ::  field_size87
REAL             , OPTIONAL ::  field87(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update88
INTEGER          , OPTIONAL ::  field_size88
REAL             , OPTIONAL ::  field88(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update89
INTEGER          , OPTIONAL ::  field_size89
REAL             , OPTIONAL ::  field89(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update90
INTEGER          , OPTIONAL ::  field_size90
REAL             , OPTIONAL ::  field90(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update91
INTEGER          , OPTIONAL ::  field_size91
REAL             , OPTIONAL ::  field91(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update92
INTEGER          , OPTIONAL ::  field_size92
REAL             , OPTIONAL ::  field92(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update93
INTEGER          , OPTIONAL ::  field_size93
REAL             , OPTIONAL ::  field93(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update94
INTEGER          , OPTIONAL ::  field_size94
REAL             , OPTIONAL ::  field94(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update95
INTEGER          , OPTIONAL ::  field_size95
REAL             , OPTIONAL ::  field95(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update96
INTEGER          , OPTIONAL ::  field_size96
REAL             , OPTIONAL ::  field96(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update97
INTEGER          , OPTIONAL ::  field_size97
REAL             , OPTIONAL ::  field97(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update98
INTEGER          , OPTIONAL ::  field_size98
REAL             , OPTIONAL ::  field98(*)
CHARACTER (len=5), OPTIONAL ::  key_to_update99
INTEGER          , OPTIONAL ::  field_size99
REAL             , OPTIONAL ::  field99(*)

IF (lhook) CALL dr_hook('INTEGRITY:UPDATE_HASH_M',zhook_in,zhook_handle)

CALL update_hash(field ,field_size ,key_to_update)
IF(PRESENT(field1)) CALL update_hash(field1,field_size1,key_to_update1)
IF(PRESENT(field2)) CALL update_hash(field2,field_size2,key_to_update2)
IF(PRESENT(field3)) CALL update_hash(field3,field_size3,key_to_update3)
IF(PRESENT(field4)) CALL update_hash(field4,field_size4,key_to_update4)
IF(PRESENT(field5)) CALL update_hash(field5,field_size5,key_to_update5)
IF(PRESENT(field6)) CALL update_hash(field6,field_size6,key_to_update6)
IF(PRESENT(field7)) CALL update_hash(field7,field_size7,key_to_update7)
IF(PRESENT(field8)) CALL update_hash(field8,field_size8,key_to_update8)
IF(PRESENT(field9)) CALL update_hash(field9,field_size9,key_to_update9)
IF(PRESENT(field10)) CALL update_hash(field10,field_size10,key_to_update10)
IF(PRESENT(field11)) CALL update_hash(field11,field_size11,key_to_update11)
IF(PRESENT(field12)) CALL update_hash(field12,field_size12,key_to_update12)
IF(PRESENT(field13)) CALL update_hash(field13,field_size13,key_to_update13)
IF(PRESENT(field14)) CALL update_hash(field14,field_size14,key_to_update14)
IF(PRESENT(field15)) CALL update_hash(field15,field_size15,key_to_update15)
IF(PRESENT(field16)) CALL update_hash(field16,field_size16,key_to_update16)
IF(PRESENT(field17)) CALL update_hash(field17,field_size17,key_to_update17)
IF(PRESENT(field18)) CALL update_hash(field18,field_size18,key_to_update18)
IF(PRESENT(field19)) CALL update_hash(field19,field_size19,key_to_update19)
IF(PRESENT(field20)) CALL update_hash(field20,field_size20,key_to_update20)
IF(PRESENT(field21)) CALL update_hash(field21,field_size21,key_to_update21)
IF(PRESENT(field22)) CALL update_hash(field22,field_size22,key_to_update22)
IF(PRESENT(field23)) CALL update_hash(field23,field_size23,key_to_update23)
IF(PRESENT(field24)) CALL update_hash(field24,field_size24,key_to_update24)
IF(PRESENT(field25)) CALL update_hash(field25,field_size25,key_to_update25)
IF(PRESENT(field26)) CALL update_hash(field26,field_size26,key_to_update26)
IF(PRESENT(field27)) CALL update_hash(field27,field_size27,key_to_update27)
IF(PRESENT(field28)) CALL update_hash(field28,field_size28,key_to_update28)
IF(PRESENT(field29)) CALL update_hash(field29,field_size29,key_to_update29)
IF(PRESENT(field30)) CALL update_hash(field30,field_size30,key_to_update30)
IF(PRESENT(field31)) CALL update_hash(field31,field_size31,key_to_update31)
IF(PRESENT(field32)) CALL update_hash(field32,field_size32,key_to_update32)
IF(PRESENT(field33)) CALL update_hash(field33,field_size33,key_to_update33)
IF(PRESENT(field34)) CALL update_hash(field34,field_size34,key_to_update34)
IF(PRESENT(field35)) CALL update_hash(field35,field_size35,key_to_update35)
IF(PRESENT(field36)) CALL update_hash(field36,field_size36,key_to_update36)
IF(PRESENT(field37)) CALL update_hash(field37,field_size37,key_to_update37)
IF(PRESENT(field38)) CALL update_hash(field38,field_size38,key_to_update38)
IF(PRESENT(field39)) CALL update_hash(field39,field_size39,key_to_update39)
IF(PRESENT(field40)) CALL update_hash(field40,field_size40,key_to_update40)
IF(PRESENT(field41)) CALL update_hash(field41,field_size41,key_to_update41)
IF(PRESENT(field42)) CALL update_hash(field42,field_size42,key_to_update42)
IF(PRESENT(field43)) CALL update_hash(field43,field_size43,key_to_update43)
IF(PRESENT(field44)) CALL update_hash(field44,field_size44,key_to_update44)
IF(PRESENT(field45)) CALL update_hash(field45,field_size45,key_to_update45)
IF(PRESENT(field46)) CALL update_hash(field46,field_size46,key_to_update46)
IF(PRESENT(field47)) CALL update_hash(field47,field_size47,key_to_update47)
IF(PRESENT(field48)) CALL update_hash(field48,field_size48,key_to_update48)
IF(PRESENT(field49)) CALL update_hash(field49,field_size49,key_to_update49)
IF(PRESENT(field50)) CALL update_hash(field50,field_size50,key_to_update50)
IF(PRESENT(field51)) CALL update_hash(field51,field_size51,key_to_update51)
IF(PRESENT(field52)) CALL update_hash(field52,field_size52,key_to_update52)
IF(PRESENT(field53)) CALL update_hash(field53,field_size53,key_to_update53)
IF(PRESENT(field54)) CALL update_hash(field54,field_size54,key_to_update54)
IF(PRESENT(field55)) CALL update_hash(field55,field_size55,key_to_update55)
IF(PRESENT(field56)) CALL update_hash(field56,field_size56,key_to_update56)
IF(PRESENT(field57)) CALL update_hash(field57,field_size57,key_to_update57)
IF(PRESENT(field58)) CALL update_hash(field58,field_size58,key_to_update58)
IF(PRESENT(field59)) CALL update_hash(field59,field_size59,key_to_update59)
IF(PRESENT(field60)) CALL update_hash(field60,field_size60,key_to_update60)
IF(PRESENT(field61)) CALL update_hash(field61,field_size61,key_to_update61)
IF(PRESENT(field62)) CALL update_hash(field62,field_size62,key_to_update62)
IF(PRESENT(field63)) CALL update_hash(field63,field_size63,key_to_update63)
IF(PRESENT(field64)) CALL update_hash(field64,field_size64,key_to_update64)
IF(PRESENT(field65)) CALL update_hash(field65,field_size65,key_to_update65)
IF(PRESENT(field66)) CALL update_hash(field66,field_size66,key_to_update66)
IF(PRESENT(field67)) CALL update_hash(field67,field_size67,key_to_update67)
IF(PRESENT(field68)) CALL update_hash(field68,field_size68,key_to_update68)
IF(PRESENT(field69)) CALL update_hash(field69,field_size69,key_to_update69)
IF(PRESENT(field70)) CALL update_hash(field70,field_size70,key_to_update70)
IF(PRESENT(field71)) CALL update_hash(field71,field_size71,key_to_update71)
IF(PRESENT(field72)) CALL update_hash(field72,field_size72,key_to_update72)
IF(PRESENT(field73)) CALL update_hash(field73,field_size73,key_to_update73)
IF(PRESENT(field74)) CALL update_hash(field74,field_size74,key_to_update74)
IF(PRESENT(field75)) CALL update_hash(field75,field_size75,key_to_update75)
IF(PRESENT(field76)) CALL update_hash(field76,field_size76,key_to_update76)
IF(PRESENT(field77)) CALL update_hash(field77,field_size77,key_to_update77)
IF(PRESENT(field78)) CALL update_hash(field78,field_size78,key_to_update78)
IF(PRESENT(field79)) CALL update_hash(field79,field_size79,key_to_update79)
IF(PRESENT(field80)) CALL update_hash(field80,field_size80,key_to_update80)
IF(PRESENT(field81)) CALL update_hash(field81,field_size81,key_to_update81)
IF(PRESENT(field82)) CALL update_hash(field82,field_size82,key_to_update82)
IF(PRESENT(field83)) CALL update_hash(field83,field_size83,key_to_update83)
IF(PRESENT(field84)) CALL update_hash(field84,field_size84,key_to_update84)
IF(PRESENT(field85)) CALL update_hash(field85,field_size85,key_to_update85)
IF(PRESENT(field86)) CALL update_hash(field86,field_size86,key_to_update86)
IF(PRESENT(field87)) CALL update_hash(field87,field_size87,key_to_update87)
IF(PRESENT(field88)) CALL update_hash(field88,field_size88,key_to_update88)
IF(PRESENT(field89)) CALL update_hash(field89,field_size89,key_to_update89)
IF(PRESENT(field90)) CALL update_hash(field90,field_size90,key_to_update90)
IF(PRESENT(field91)) CALL update_hash(field91,field_size91,key_to_update91)
IF(PRESENT(field92)) CALL update_hash(field92,field_size92,key_to_update92)
IF(PRESENT(field93)) CALL update_hash(field93,field_size93,key_to_update93)
IF(PRESENT(field94)) CALL update_hash(field94,field_size94,key_to_update94)
IF(PRESENT(field95)) CALL update_hash(field95,field_size95,key_to_update95)
IF(PRESENT(field96)) CALL update_hash(field96,field_size96,key_to_update96)
IF(PRESENT(field97)) CALL update_hash(field97,field_size97,key_to_update97)
IF(PRESENT(field98)) CALL update_hash(field98,field_size98,key_to_update98)
IF(PRESENT(field99)) CALL update_hash(field99,field_size99,key_to_update99)

IF (lhook) CALL dr_hook('INTEGRITY:UPDATE_HASH_M',zhook_out,zhook_handle)

END SUBROUTINE


SUBROUTINE update_hash(field,field_size,key_to_update)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE proc_info_mod,     ONLY : me,n_proc

!#if defined (1) && (EG_INTEGRITY_TEST)
!USE IFCORE
!#endif

USE timestep_mod, ONLY : timestep_number

IMPLICIT NONE
CHARACTER (len=9) filename
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (len=5)  key_to_update
INTEGER field_size,ierr
REAL    field(*),field_max,field_min
CHARACTER(len=65) hash_to_update

INTEGER entry
!
!  Purpose: Defines unit numbers relevant to history file
!           and variables used to hold the logical to physical
!           file associations made within the model
!
!  Logical Filenames used in the model
!
      CHARACTER(LEN=256) hkfile,ppxref,config,stashctl,namelist,output,      &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,ftxx,    &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar
     

!
      CHARACTER(LEN=256) MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                                 ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        RSUB_UNIT,                                                &
                                 ! File indicating whether resub required
     &        XHIST_UNIT,                                               &
                                 ! Main history file unit
     &        THIST_UNIT,                                               &
                                 ! Backup history file unit
     &        HKFILE_UNIT,                                              &
                                 ! Operational houskeeping file unit    
     &        EG_UNIT            ! ENDGame diagnostics/info unit
!
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(RSUB_UNIT =10)
      PARAMETER(XHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)

!
! Parameters specifying unit numbers relevant to ENDGame diagnostics
!
      PARAMETER(EG_UNIT  = 55)

! UKCA unit numbers

      INTEGER, PARAMETER :: ukcafjxx_unit=170 ! Fast-J(X) cross section data
      INTEGER, PARAMETER :: ukcafjsc_unit=171 ! Fast-JX scattering data
      INTEGER, PARAMETER :: ukca2do3_unit=172 ! 2D top boundary O3 data 
      INTEGER, PARAMETER :: ukca2ch4_unit=173 ! 2D top boundary CH4 data
      INTEGER, PARAMETER :: ukca2noy_unit=174 ! 2D top boundary NOY data
      INTEGER, PARAMETER :: ukca2pho_unit=175 ! 2D photolysis input data
      INTEGER, PARAMETER :: ukcastrd_unit=176 ! Stratospheric model radiation field. 
      INTEGER, PARAMETER :: ukcasto3_unit=177 ! Strat standard atmosphere T and O3.
      INTEGER, PARAMETER :: ukcastar_unit=178 ! Stratospheric sulfate aerosol climatology 
      INTEGER, PARAMETER :: ukcafjar_unit=179 ! Sulfate aerosol cliamtology for Fast-JX
! Text output file for STASH-related information is assigned to UNIT 200

!
! Namelist of all permissible logical files.
!
      NAMELIST / nlcfiles /                                             &
                   hkfile,ppxref,config,stashctl,namelist,output,       &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,         &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(RSUB       ,MODEL_FT_UNIT(10) ), &
     &(XHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(ICECALVE  ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &                                (ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(GENLAND   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
                                      (lwspectd   ,model_ft_unit(80) ), &
      (surgeou1  ,model_ft_unit(81) ),(surgeout   ,model_ft_unit(82) ), &
      (ppscreen  ,model_ft_unit(83) ),(ppsmc      ,model_ft_unit(84) ), &
      (wfout     ,model_ft_unit(85) ),(uarsout1   ,model_ft_unit(86) ), &
      (uarsout2  ,model_ft_unit(87) ),(icefout    ,model_ft_unit(88) ), &
      (mosout    ,model_ft_unit(89) ),(vert_lev   ,model_ft_unit(90) ), &
      (sstout    ,model_ft_unit(91) ),(siceout    ,model_ft_unit(92) ), &
      (curntout  ,model_ft_unit(93) ),(flxcrout   ,model_ft_unit(94) ), &
      (dmsconc   ,model_ft_unit(95) ),(orog       ,model_ft_unit(96) ), &
      (transp    ,model_ft_unit(97) ),(olabcin    ,model_ft_unit(98) ), &
      (ocndepth  ,model_ft_unit(99) ),                                  &
      (foamout1  ,model_ft_unit(100)),(foamout2   ,model_ft_unit(101)), &
      (cxbkgerr  ,model_ft_unit(102)),(rfmout     ,model_ft_unit(103)), &
      (idealise  ,model_ft_unit(106)),(tdf_dump   ,model_ft_unit(107)), &
      (iau_inc   ,model_ft_unit(108)),(murkfile   ,model_ft_unit(109)), &
      (sulpemis  ,model_ft_unit(110)),(usrancil   ,model_ft_unit(111)), &
      (usrmulti  ,model_ft_unit(112)),(ousrancl   ,model_ft_unit(113)), &
      (ousrmult  ,model_ft_unit(114)),(so2natem   ,model_ft_unit(115)), &
      (chemoxid  ,model_ft_unit(116)),(aerofcg    ,model_ft_unit(117)), &
      (co2emits  ,model_ft_unit(118)),(tppsozon   ,model_ft_unit(119)), &
      (landfrac  ,model_ft_unit(120)),(wlabcou1   ,model_ft_unit(121)), &
      (wlabcou2  ,model_ft_unit(122)),(wlabcou3   ,model_ft_unit(123)), &
      (wlabcou4  ,model_ft_unit(124)),(alabcin1   ,model_ft_unit(125)), &
      (alabcin2  ,model_ft_unit(126)),                                  &
      (ocffemis  ,model_ft_unit(128)),(horzgrid   ,model_ft_unit(129)), &
      (surfemis  ,model_ft_unit(130)),(aircrems   ,model_ft_unit(131)), &
      (stratems  ,model_ft_unit(132)),(extraems   ,model_ft_unit(133)), &
      (radonems  ,model_ft_unit(134)),(fracinit   ,model_ft_unit(135)), &
      (veginit   ,model_ft_unit(136)),(disturb    ,model_ft_unit(137)), &
      (cached    ,model_ft_unit(138)),(sootemis   ,model_ft_unit(139)), &
      (alabcou1  ,model_ft_unit(140)),(alabcou2   ,model_ft_unit(141)), &
      (alabcou3  ,model_ft_unit(142)),(alabcou4   ,model_ft_unit(143)), &
      (alabcou5  ,model_ft_unit(144)),(alabcou6   ,model_ft_unit(145)), &
      (alabcou7  ,model_ft_unit(146)),(alabcou8   ,model_ft_unit(147)), &
      (cariolo3  ,model_ft_unit(148)),(rpseed     ,model_ft_unit(149)), &
      (ppvar     ,model_ft_unit(150)),(pp10       ,model_ft_unit(151)), &
      (icfile    ,model_ft_unit(152)),(var_grid   ,model_ft_unit(153)), &
      (arclbiog  ,model_ft_unit(154)),(arclbiom   ,model_ft_unit(155)), &
      (arclblck  ,model_ft_unit(156)),(arclsslt   ,model_ft_unit(157)), &
      (arclsulp  ,model_ft_unit(158)),(arcldust   ,model_ft_unit(159)), &
      (arclocff  ,model_ft_unit(160)),(arcldlta   ,model_ft_unit(161)), &
      (topmean   ,model_ft_unit(162)),(topstdev   ,model_ft_unit(163)), &
      (ppmbc     ,model_ft_unit(164)),(ukcaprec   ,model_ft_unit(165)), &
      (ukcaacsw  ,model_ft_unit(166)),(ukcaaclw   ,model_ft_unit(167)), &
      (ukcacrsw  ,model_ft_unit(168)),(ukcacrlw   ,model_ft_unit(169)), &
      (ukcafjxx  ,model_ft_unit(170)),(ukcafjsc   ,model_ft_unit(171)), &
      (ukca2do3  ,model_ft_unit(172)),(ukca2ch4   ,model_ft_unit(173)), &
      (ukca2noy  ,model_ft_unit(174)),(ukca2pho   ,model_ft_unit(175)), &
      (ukcastrd  ,model_ft_unit(176)),(ukcasto3   ,model_ft_unit(177)), &
      (ukcastar  ,model_ft_unit(178)),(ukcafjar   ,model_ft_unit(179))
! Text output file for STASH-related information is assigned to UNIT 200


!#if defined (IBM) && (EG_INTEGRITY_TEST)
!        include 'fexcp.h'
!#endif

IF (lhook) CALL dr_hook('INTEGRITY:UPDATE_HASH',zhook_in,zhook_handle)

hash_to_update="                                                                  "

entry=1
DO WHILE (key(entry) /= key_to_update .AND. entry.lt.registry_length)
  entry=entry+1
END DO

IF(entry.le.registry_length .and. key(entry) == key_to_update) THEN




  IF(ANY(field(1:field_size)/=field(1:field_size))) THEN
  

!#if defined (IBM) && (EG_INTEGRITY_TEST)
!        call xl__trbk
!#endif
!#if defined (1) && (EG_INTEGRITY_TEST)
!CALL TRACEBACKQQ(user_exit_code=0)
!#endif

!    cmessage = ' NaN detected in update of ' // TRIM(key_to_update)

!    CALL ereport(routinename,1,cmessage)

  ENDIF

  field_max=maxval(field(1:field_size))
  field_min=minval(field(1:field_size))

  CALL gc_rmax(1,n_proc,ierr,field_max)
  CALL gc_rmin(1,n_proc,ierr,field_min)

  ! DEPENDS ON: eg_hash
  CALL eg_hash(field,8*field_size,hash_to_update)

  hash(entry)=hash_to_update


    IF (me == 0) WRITE (av_out,fmt='(A,I4.4,2A,I12,2E16.8,2A)') 'updating PE '&
         ,me,&
        ': ',key(entry),field_size,field_min,   &
       field_max,' ', hash_to_update(1:64)


  ! DEPENDS ON: eg_hash
  CALL eg_hash(hash,65*(last_entry),hash(registry_length+1))

  IF (integrity_file) THEN
    WRITE (filename,fmt='(I4.4,A,I4.4)') me,'-',timestep_number
    OPEN(eg_unit,file=filename, position='APPEND')
    WRITE(eg_unit,'(A,A)') key(entry),hash_to_update(1:63)
    CLOSE(eg_unit)
  END IF
ELSE
! not found for update, so add it as new
  CALL add_hash(field,field_size,key_to_update)
END IF

IF (lhook) CALL dr_hook('INTEGRITY:UPDATE_HASH',zhook_out,zhook_handle)

END SUBROUTINE


SUBROUTINE list_hash()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE proc_info_mod,     ONLY : me

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER i
CHARACTER(len=65) hashvalue
IF (lhook) CALL dr_hook('INTEGRITY:LIST_HASH',zhook_in,zhook_handle)

IF (me==0) THEN
  WRITE(av_out,fmt='(A)') &
      '========================================================================'
  WRITE(av_out,fmt='(A)') &
      ' ENDGame data integrity test - checksum table                           '
  WRITE(av_out,fmt='(A)') &
      '                                                                        '
  WRITE(av_out,fmt='(A)') &
      ' key:   checksum:                                                       '
  WRITE(av_out,fmt='(A)') &
      ' -----------------------------------------------------------------------'

  DO i=1,last_entry
    hashvalue=hash(i)
    WRITE(av_out,fmt='(4A)') ' ',key(i),'  ',hashvalue(1:64)
  END DO
  WRITE(av_out,fmt='(A)') &
      '========================================================================'
END IF

IF (lhook) CALL dr_hook('INTEGRITY:LIST_HASH',zhook_out,zhook_handle)

END SUBROUTINE

SUBROUTINE global_hash(field,key,dim_in,fld_type,unit_in)

USE atm_fields_bounds_mod, ONLY : array_dims
USE proc_info_mod,         ONLY : me,n_proc,n_procx,n_procy,global_row_length

IMPLICIT NONE

TYPE (array_dims) :: dim_in

REAL ::       field(dim_in%i_start:dim_in%i_end, &
                    dim_in%j_start:dim_in%j_end, &
                    dim_in%k_start:dim_in%k_end)

CHARACTER (len=10) :: key

INTEGER            :: fld_type        ! type (p,u or v) of grid              
INTEGER, OPTIONAL  :: unit_in

IF (PRESENT(unit_in)) global_unt = unit_in

IF(n_proc>512.or.global_row_length>640) THEN
!
! liklely to be a too big job for frequent all_reduce of a 3D field
! as required for this function. Therefore switching to a cheaper version
! if more than 1024 cpus or job > N320

! This should possibly be turned into a runtime switch at some point!

  CALL global_hash1(field,key,dim_in,fld_type)

ELSE

  CALL global_hash2(field,key,dim_in,fld_type)

END IF

END SUBROUTINE

SUBROUTINE global_hash1(field,key,dim_in,fld_type,unit_in)

!
! This version does compute a hash on each cpu and then reports
! a hash of hashes. This of course only compares when the cpu
! count and decomposition is not changed, but is substantially faster
! and less memory hungy for large jobs
!

USE ereport_mod,           ONLY : ereport
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE UM_ParVars
USE atm_fields_bounds_mod, ONLY : array_dims
USE proc_info_mod,         ONLY : me,n_proc,n_procx,n_procy
USE MPL

IMPLICIT NONE

! ARGUMENTS

TYPE (array_dims)      :: dim_in

REAL                   :: field(dim_in%i_start:dim_in%i_end, &
                                dim_in%j_start:dim_in%j_end, &
                                dim_in%k_start:dim_in%k_end)
CHARACTER (len=10)     :: key

INTEGER                :: fld_type      ! type (p,u or v) of grid  
INTEGER, OPTIONAL      :: unit_in


! LOCAL

INTEGER                :: local_row_len  
INTEGER                :: local_rows     
INTEGER                :: global_row_len        
INTEGER                :: global_rows          
INTEGER                :: icode            

CHARACTER (len=256)    :: CMessage       

CHARACTER (len=65)     :: ghash(n_proc),pe_hash

INTEGER                :: k
INTEGER                :: ier,array_size,root
INTEGER                :: gcom_mpi_comm_world,istat

!xternal mpi_gather

IF (PRESENT(unit_in)) global_unt = unit_in

CALL eg_hash(field,8*SIZE(field),pe_hash)

CALL gc_get_communicator(gcom_mpi_comm_world, istat)

root = 0
array_size = 65

CALL MPL_GATHER ( pe_hash,             &
                  array_size,          &
                  MPL_CHARACTER,       &
                  ghash,               &
                  array_size,          &
                  MPL_CHARACTER,       &
                  root,                &
                  gcom_MPI_COMM_WORLD, &
                  ier)


IF (me == 0) THEN
! DEPENDS ON: eg_hash
  CALL eg_hash(ghash,65*n_proc,pe_hash)

  WRITE(global_unt,fmt='(3A)') 'gHash of (only valid for this decomp)'&
                              ,key,pe_hash
END IF

END SUBROUTINE


SUBROUTINE global_hash2(field,key,dim_in,fld_type,unit_in)

!
! Reconstruct the global field on CPU 0, then compute the hash. This
! will be the same, irrespective of cpu count and decomposition, 
! enabeling bitcomparision tests
!

USE ereport_mod,           ONLY : ereport
USE UM_ParVars
USE atm_fields_bounds_mod, ONLY : array_dims
USE proc_info_mod,         ONLY : me,n_proc,n_procx,n_procy


IMPLICIT NONE


! ARGUMENTS
TYPE (array_dims)  :: dim_in

REAL               :: field(dim_in%i_start:dim_in%i_end, &
                            dim_in%j_start:dim_in%j_end, &
                            dim_in%k_start:dim_in%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (len=10) :: key
INTEGER            :: fld_type        ! type (p,u or v) of grid        
INTEGER, OPTIONAL  :: unit_in

! LOCAL VARIABLES

INTEGER            :: local_row_len   ! length of rows in
                                      ! local part of field
INTEGER            :: local_rows      ! number of rows in
                                      ! local part of field
INTEGER            :: global_row_len  ! length of rows in
                                      ! global field       
INTEGER            :: global_rows     ! number of rows in
                                      ! global field       
      
INTEGER            :: icode           ! out return code   

REAL, ALLOCATABLE  :: local_field(:,:)
REAL, ALLOCATABLE  :: global_field(:,:)

CHARACTER (len=256):: CMessage       ! Error message if return code >0

CHARACTER (len=65) :: ghash(dim_in%k_start:dim_in%k_end+1)

INTEGER            :: k
INTEGER            :: unt

IF (PRESENT(unit_in)) global_unt = unit_in


IF (lhook) CALL dr_hook('INTEGRITY:GLOBAL_HASH',zhook_in,zhook_handle)

unt = 0
IF (PRESENT(unit_in)) unt = unit_in

local_row_len = dim_in%i_end-dim_in%i_start+1-2*dim_in%halo_i
local_rows    = dim_in%j_end-dim_in%j_start+1-2*dim_in%halo_j

ALLOCATE (local_field(local_row_len,local_rows))

global_row_len = local_row_len
CALL gc_isum (1,n_proc,icode,  global_row_len)
global_rows    = local_rows   
CALL gc_isum (1,n_proc,icode,  global_rows   )

global_row_len = global_row_len / n_procy
global_rows    = global_rows    / n_procx

IF (me == 0) THEN
  ALLOCATE (global_field(global_row_len,global_rows))
ELSE
  ALLOCATE (global_field(1,1))
END IF

DO k=dim_in%k_start,dim_in%k_end

  local_field(:,:) = field(dim_in%i_start+dim_in%halo_i:              &
                           dim_in%i_end  -dim_in%halo_i               &
                          ,dim_in%j_start+dim_in%halo_j:              &
                           dim_in%j_end  -dim_in%halo_j               &
                          ,k)

!DEPENDS ON:gather_field
  CALL  gather_field(local_field,    global_field,                    &
                        local_row_len,  local_rows,                   &
                        global_row_len, global_rows,                  &
                        fld_type,      halo_type_no_halo,             &
                        0,     GC_ALL_PROC_GROUP,                     &
                        icode,          cmessage)          

  IF (me == 0)                                                        &
  ! DEPENDS ON: eg_hash
    CALL eg_hash(global_field,8*SIZE(global_field),ghash(k))

END DO

IF (me == 0) THEN

! DEPENDS ON: eg_hash
  CALL eg_hash(ghash,65*(dim_in%k_end-dim_in%k_start+1),              &
                     ghash(dim_in%k_end+1))

  WRITE(global_unt,fmt='(3A)') 'gHash of ',key,ghash(dim_in%k_end+1)

END IF

DEALLOCATE(global_field)
DEALLOCATE(local_field)

IF (lhook) CALL dr_hook('INTEGRITY:GLOBAL_HASH',zhook_out,zhook_handle)

END SUBROUTINE

END MODULE
