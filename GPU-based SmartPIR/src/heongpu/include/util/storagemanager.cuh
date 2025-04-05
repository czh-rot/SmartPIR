// Copyright 2024 Alişah Özcan
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Alişah Özcan

#ifndef STORAGE_MANAGER_H
#define STORAGE_MANAGER_H

#include "common.cuh"
#include "nttparameters.cuh"
#include <iostream>
#include <stdexcept>
#include <vector>

namespace heongpu
{
    /**
     * @enum storage_type
     * @brief Specifies the storage type for objects.
     *
     * - `HOST`: Indicates that the object is stored in the host memory.
     * - `DEVICE`: Indicates that the object is stored in the device memory.
     */
    enum class storage_type : std::uint8_t
    {
        HOST = 0x1, // Host memory storage
        DEVICE = 0x2 // Device memory storage
    };

    /**
     * @struct ExecutionOptions
     * @brief Manages execution options, including CUDA stream configuration,
     * storage type selection, and the behavior regarding initial data location.
     */
    struct ExecutionOptions
    {
        cudaStream_t stream_ =
            cudaStreamDefault; // CUDA stream to be used for execution. Defaults
                               // to cudaStreamDefault.
        storage_type storage_ =
            storage_type::DEVICE; // Default storage type for the operation.
                                  // Defaults to DEVICE.
        bool keep_initial_condition_ =
            true; // Determines whether to maintain the initial data location.
                  // Defaults to true.

        /**
         * @brief Sets the CUDA stream to be used for execution.
         *
         * This method allows you to specify the CUDA stream that will be used
         * during the execution of operations. By default, the stream is set to
         * `cudaStreamDefault`.
         *
         * @param stream The CUDA stream to set.
         * @return ExecutionOptions& A reference to the updated ExecutionOptions
         * object, allowing method chaining.
         */
        ExecutionOptions& set_stream(cudaStream_t stream)
        {
            stream_ = stream;
            return *this;
        }

        /**
         * @brief Configures the storage type for execution.
         *
         * Sets the storage type where the operation's data will be stored. This
         * can typically be either `storage_type::DEVICE` or
         * `storage_type::HOST`.
         *
         * @param storage The storage type to set (e.g., DEVICE or HOST).
         * @return ExecutionOptions& A reference to the updated ExecutionOptions
         * object, allowing method chaining.
         */
        ExecutionOptions& set_storage_type(storage_type storage)
        {
            storage_ = storage;
            return *this;
        }

        /**
         * @brief Specifies whether to preserve the initial data location.
         *
         * Controls whether the initial data location (e.g., host or device) is
         * preserved during execution. If set to `true`, the data will be
         * returned to its initial location after the operation.
         *
         * @param keep_initial_condition Boolean flag to determine whether to
         * keep the initial location.
         * @return ExecutionOptions& A reference to the updated ExecutionOptions
         * object, allowing method chaining.
         */
        ExecutionOptions& set_initial_location(bool keep_initial_condition)
        {
            keep_initial_condition_ = keep_initial_condition;
            return *this;
        }
    };

    /**
     * @brief Manages the input storage and conditionally transfers data to the
     * appropriate location (e.g., device or host) before executing a function
     * on the object.
     *
     * @tparam T The type of the object to manage.
     * @tparam F The callable type of the function to execute on the object.
     * @param object The object whose storage is being managed.
     * @param function The callable function that will operate on the object.
     * @param options A set of execution options defining storage behavior and
     * stream configuration.
     * @param is_input_output_same A flag indicating whether the input and
     * output storage should remain consistent.
     */
    template <typename T, typename F>
    void input_storage_manager(T& object, F function, ExecutionOptions options,
                               bool is_input_output_same)
    {
        storage_type initial_condition = object.storage_type_;

        if (!object.is_on_device())
        {
            if (options.keep_initial_condition_)
            {
                object.copy_to_device(options.stream_);
            }
            else
            {
                object.store_in_device(options.stream_);
            }
        }

        function(object);

        if (is_input_output_same)
        {
            if (options.keep_initial_condition_ &&
                (initial_condition == storage_type::HOST) &&
                (options.storage_ == storage_type::DEVICE))
            {
                object.remove_from_host();
            }
        }
        else
        {
            if (options.keep_initial_condition_)
            {
                if (initial_condition == storage_type::HOST)
                {
                    object.remove_from_device(options.stream_);
                }
            }
            else
            {
                if (options.storage_ == storage_type::DEVICE)
                {
                    object.store_in_device(options.stream_);
                }
                else if (options.storage_ == storage_type::HOST)
                {
                    object.store_in_host(options.stream_);
                }
                else
                {
                    throw std::invalid_argument("Invalid storage type!");
                }
            }
        }
    }

    /**
     * @brief Recursively manages the input storage for a vector of objects and
     * conditionally transfers data to the appropriate location (e.g., device or
     * host) before executing a function on each object.
     *
     * @tparam T The type of the objects in the vector to manage.
     * @tparam F The callable type of the function to execute on each object.
     * @param objects A vector of objects whose storage is being managed.
     * @param function The callable function that will operate on each object in
     * the vector.
     * @param options A set of execution options defining storage behavior and
     * stream configuration.
     * @param is_input_output_same A flag indicating whether the input and
     * output storage should remain consistent for each object.
     */
    template <typename T, typename F>
    void input_vector_storage_manager(std::vector<T>& objects, F function,
                                      ExecutionOptions options,
                                      bool is_input_output_same)
    {
        std::vector<storage_type> initial_conditions(objects.size());

        for (int i = 0; i < objects.size(); i++)
        {
            if (!objects[i].is_on_device())
            {
                if (options.keep_initial_condition_)
                {
                    objects[i].copy_to_device(options.stream_);
                }
                else
                {
                    objects[i].store_in_device(options.stream_);
                }
            }
        }

        function(objects);

        for (int i = 0; i < objects.size(); i++)
        {
            if (is_input_output_same)
            {
                if (options.keep_initial_condition_ &&
                    (initial_conditions[i] == storage_type::HOST) &&
                    (options.storage_ == storage_type::DEVICE))
                {
                    objects[i].remove_from_host();
                }
            }
            else
            {
                if (options.keep_initial_condition_)
                {
                    if (initial_conditions[i] == storage_type::HOST)
                    {
                        objects[i].remove_from_device(options.stream_);
                    }
                }
                else
                {
                    if (options.storage_ == storage_type::DEVICE)
                    {
                        objects[i].store_in_device(options.stream_);
                    }
                    else if (options.storage_ == storage_type::HOST)
                    {
                        objects[i].store_in_host(options.stream_);
                    }
                    else
                    {
                        throw std::invalid_argument("Invalid storage type!");
                    }
                }
            }
        }
    }

    /**
     * @brief Manages the output storage of an object after executing a
     * function.
     *
     * @tparam T The type of the object to manage.
     * @tparam F The callable type of the function to execute on the object.
     * @param object The object whose output storage is being managed.
     * @param function The callable function that will operate on the object.
     * @param options A set of execution options defining storage behavior and
     * stream configuration.
     */
    template <typename T, typename F>
    void output_storage_manager(T& object, F function, ExecutionOptions options)
    {
        function(object);

        if (options.storage_ == storage_type::DEVICE)
        {
            object.store_in_device(options.stream_);
        }
        else if (options.storage_ == storage_type::HOST)
        {
            object.store_in_host(options.stream_);
        }
        else
        {
            throw std::invalid_argument("Invalid storage type!");
        }
    }

} // namespace heongpu
#endif // STORAGE_MANAGER_H