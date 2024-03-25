// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

void htrun(void (*lcallfunc)(int, void *), void *luserdata, int v0, int v1, bool should_parallelise);

#endif // THREAD_POOL_HPP
