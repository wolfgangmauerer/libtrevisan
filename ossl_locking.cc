/* This file is part of libtrevisan, a modular implementation of
   Trevisan's randomness extraction construction.

   Copyright (C) 2011-2012, Wolfgang Mauerer <wm@linux-kernel.net>

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with libtrevisan. If not, see <http://www.gnu.org/licenses/>. */

#include<iostream>
#include<vector>
#include<mutex>
#include<openssl/crypto.h>
#include <pthread.h>

struct CRYPTO_dynlock_value {
	std::mutex mutex;
};

using namespace std;
static vector<class std::mutex*> ossl_mutexes;

#if OPENSSL_VERSION_NUMBER >= 0x01000000f
static void ossl_thread_id(CRYPTO_THREADID *id) {
#if EXPENSIVE_SANITY_CHECKS
	cerr << "pthreads_thread_id=" << pthread_self() << " determined" << endl;
#endif
	CRYPTO_THREADID_set_pointer(id, (void*) pthread_self());
}
#else
static unsigned long ossl_thread_id() {
#if EXPENSIVE_SANITY_CHECKS
	cerr << "pthreads_thread_id=" << pthread_self() << " determined" << endl;
#endif
	return (unsigned long)pthread_self();
}
#endif

static void ossl_locking_callback(int mode, int type, const char *file, int line) {
#if EXPENSIVE_SANITY_CHECKS
	cerr << "thread=" << pthread_self() << ", " << "mode="
	     << ((mode&CRYPTO_LOCK)?"l":"u") << ", lock="
	     << ((type&CRYPTO_READ)?"r":"w") << ", file=" << file
	     << ", line=" << line << endl;
#endif
	if (mode & CRYPTO_LOCK) {
		ossl_mutexes[type]->lock();
	} else {
		ossl_mutexes[type]->unlock();
	}
}

static void ossl_dynlock_lock(int mode, struct CRYPTO_dynlock_value *lock,
			      const char *file, int line) {
	cout << "OPENSSL: dynlock operation @" << lock << endl;
	if (mode & CRYPTO_LOCK) {
		lock->mutex.lock();
	} else {
		lock->mutex.unlock();
	}
}

static struct CRYPTO_dynlock_value *ossl_dynlock_create(const char *file, int line) {
	struct CRYPTO_dynlock_value *lock = new CRYPTO_dynlock_value;
	cout << "OPENSSL: Created new dynamic lock @" << lock << endl;

	return lock;
}

static void ossl_dynlock_destroy(struct CRYPTO_dynlock_value *lock, const char *file,
				int line) {
	delete lock;
}

void init_openssl_locking() {
	unsigned num_locks = CRYPTO_num_locks();

	ossl_mutexes.reserve(num_locks);
	for (unsigned i=0; i < num_locks; i++) {
		ossl_mutexes.push_back(new std::mutex());
	}

#if OPENSSL_VERSION_NUMBER >= 0x01000000f
	CRYPTO_THREADID_set_callback(ossl_thread_id);
#else
	CRYPTO_set_id_callback(ossl_thread_id);
#endif
	CRYPTO_set_locking_callback(ossl_locking_callback);
	CRYPTO_set_dynlock_create_callback(ossl_dynlock_create);
	CRYPTO_set_dynlock_destroy_callback(ossl_dynlock_destroy);
	CRYPTO_set_dynlock_lock_callback(ossl_dynlock_lock);
}
