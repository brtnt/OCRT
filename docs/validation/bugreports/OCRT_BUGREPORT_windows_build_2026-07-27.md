# OCRT 버그 보고 — 윈도우(MinGW) 빌드가 링크 단계에서 실패한다

- 작성일: 2026-07-27
- 대상 판본: `OCRT-v1.2-2026-07-26-KST-stage2-eap-phyto-scattering-disabled-integrated`
- 앞선 보고서(엽록소 파장 상한 850 nm)와는 별개 사안이다.

---

## 1. 증상

윈도우용 컴파일러로 빌드하면 컴파일은 끝나지만 링크에서 실패한다.

```
src/rt_solver.c: warning: implicit declaration of function 'setenv'
src/rt_solver.c: warning: implicit declaration of function 'unsetenv'
x86_64-w64-mingw32-ld: rt_solver.o: undefined reference to `setenv'
x86_64-w64-mingw32-ld: rt_solver.o: undefined reference to `unsetenv'
collect2: error: ld returned 1 exit status
```

빠지는 것은 `setenv` 와 `unsetenv` 두 함수뿐이다. 나머지 34개 소스는 모두 정상 컴파일된다.

## 2. 원인

`src/main.c:15~20` 에 이미 윈도우용 시늉 함수가 있다.

```c
#ifdef _WIN32
/* Windows (MinGW) portability shims: POSIX setenv/mkdir differ on Windows. Linux build unaffected. */
#include <direct.h>
#define mkdir(path, mode) _mkdir(path)
static int setenv(const char *name, const char *value, int overwrite) { (void)overwrite; return _putenv_s(name, value); }
#endif
```

즉 **윈도우 빌드는 원래 지원 대상이다.** 다만 두 가지가 빠져 있다.

1. 이 시늉 함수가 `static` 이라 `main.c` 안에서만 보인다. `src/rt_solver.c` 도 `setenv` 를 쓰는데 그 파일에는 정의가 없다.
2. `unsetenv` 의 시늉 함수는 아예 없다. `src/rt_solver.c` 가 이 함수를 쓴다.

`setenv`/`unsetenv` 사용 지점은 다음과 같다.

| 파일 | 줄 | 함수 |
|---|---|---|
| `src/main.c` | 1157, 1915 | setenv (자체 시늉 함수로 해결됨) |
| `src/rt_solver.c` | 4612, 5129, 5325, 5909 | setenv (해결 안 됨) |
| `src/rt_solver.c` | 4788, 5893, 5910 | unsetenv (해결 안 됨) |

## 3. 영향

- 윈도우에서 이 판본을 빌드할 수 없다. 리눅스나 WSL 없이는 실행이 불가능하다.
- `-Wimplicit-function-declaration` 경고가 이미 나므로, 정순(STRICT) 빌드 방침과도 어긋난다. 리눅스에서는 `setenv` 가 표준 헤더에 있어 경고가 나지 않아 그동안 드러나지 않았다.

## 4. 제안하는 조치

1. 시늉 함수를 공용 헤더로 옮긴다. 예를 들어 `src/rt_win_compat.h` 를 만들어 `#ifdef _WIN32` 안에서 `static inline` 으로 `setenv` 와 `unsetenv` 를 정의하고, 이를 쓰는 모든 파일에서 포함한다.
2. `main.c` 안의 개별 시늉 함수는 지우고 공용 헤더를 쓰도록 바꾼다. 정의가 두 군데 있으면 나중에 또 어긋난다.
3. `unsetenv` 는 윈도우에서 `_putenv_s(name, "")` 로 구현하면 된다. 빈 문자열을 넣으면 변수가 지워진다.
4. 윈도우 빌드를 정기 점검 항목에 넣는다. 교차 컴파일러로 링크까지만 확인해도 이런 종류의 결함은 걸러진다.

## 5. 임시 조치 기록

이 보고서와 함께 전달한 `ocrt_win_compat.c` 는 기존 소스를 고치지 않고 두 함수만 채워 넣는 임시 파일이다. 이 파일을 빌드 목록에 더하면 링크가 통과한다. `main.c` 의 `static` 정의는 그 파일 안에서만 유효하므로 충돌하지 않는다.

이 임시 파일로 만든 윈도우 실행 파일의 확인 결과는 다음과 같다.

- 필요한 DLL 은 `KERNEL32.dll` 과 `msvcrt.dll` 뿐이다. 정적 링크(`-static`)로 만들어 별도 런타임 배포가 필요 없다.
- **다만 실제 윈도우에서 실행해 보지는 못했다.** 만든 환경에 윈도우 실행 수단이 없었다. 사용 전에 기준 사례 대조가 필요하다.

## 6. 참고 — 아키텍처 플래그 확인

임시 윈도우 빌드는 `-march=x86-64-v3` 를 썼다. 배포 기본값인 `-march=cascadelake` 는 AVX-512 를 요구해 그 명령어가 없는 CPU 에서 실행 자체가 실패하기 때문이다.

두 플래그가 결과를 바꾸는지 리눅스에서 확인했다. 아래 두 기준 사례에서 `cascadelake` 빌드와 `x86-64-v3` 빌드의 출력이 **비트 단위로 같다.**

| 사례 | I | Q | U |
|---|---|---|---|
| 대기만, 흑색 프레넬, 555 nm, 태양천정각 25도, 관측천정각 30도, 방위각 90도, 풍속 5, 광학두께 0.1, r50f05v01 | 4.4991861596e-02 | −9.0983939253e-04 | −7.5031808484e-03 |
| 위와 같은 기하에 해수 결합(엽록소 1.0, 총부유물 1.5, CDOM 0.035) | 1.2313717867e-01 | −6.3039574299e-04 | −1.0877124534e-02 |

결합 사례의 산란차수는 131 이고 수렴 표시는 1 이다.
