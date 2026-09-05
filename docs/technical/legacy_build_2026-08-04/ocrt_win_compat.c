/* 윈도우(MinGW) 빌드용 임시 호환 파일.
   src/main.c 에는 이미 static setenv 시늉 함수가 있으나 파일 안에서만 보이고,
   src/rt_solver.c 가 쓰는 setenv/unsetenv 는 정의가 없어 링크가 실패한다.
   이 파일은 기존 소스를 고치지 않고 그 두 함수만 채운다.
   정식 수정은 OCRT 개발 세션에서 공용 헤더로 옮기는 것이 맞다. */
#ifdef _WIN32
#include <stdlib.h>
int setenv(const char *name, const char *value, int overwrite)
{
    (void)overwrite;
    return _putenv_s(name, value);
}
int unsetenv(const char *name)
{
    return _putenv_s(name, "");
}
#endif
