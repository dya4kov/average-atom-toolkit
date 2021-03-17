set(BESSEL_SRC "" PARENT_SCOPE)
install(FILES I0.h I1.h In.h Inu.h
        J0.h J1.h Jn.h Jnu.h
        K0.h K1.h Kn.h Knu.h
        Y0.h Y1.h Yn.h Ynu.h
        DESTINATION include/${PROJECT_NAME}/${DIR_NAME}/bessel
)
