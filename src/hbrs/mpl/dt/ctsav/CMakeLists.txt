# Copyright (c) 2019 Jakob Meng, <jakobmeng@web.de>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

#################### build ####################

target_sources(hbrs_mpl PRIVATE
    impl.cpp)

#################### tests ####################

hbrs_mpl_add_test(dt_ctsav "test.cpp")
