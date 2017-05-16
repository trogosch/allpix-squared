/**
 * @file
 * @brief Base of detector models
 *
 * @copyright MIT License
 */

/**
 * @defgroup DetectorModels Detector models
 * @brief Collection of global detector models supported by the framework
 */

#ifndef ALLPIX_GEOMETRY_DETECTOR_MODEL_H
#define ALLPIX_GEOMETRY_DETECTOR_MODEL_H

#include <string>
#include <utility>

#include <Math/Point3D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>

namespace allpix {
    /**
     * @ingroup DetectorModels
     * @brief Base of all detector models
     *
     * Implements the minimum required for a detector model. A model always has a sensitive device (referred as the sensor).
     * This sensor has the following properties:
     * - A size in three dimensions
     * - A set of two dimensional replicated blocks called pixels
     * - An offset of the sensor minimum from the center of the first pixel (the origin in the local coordinate system)
     * - The coordinates of the position that is used to place the detector in the global frame (the rotation center)
     */
    class DetectorModel {
    public:
        /**
         * @brief Constructs a detector model of a certain type
         * @param type Unique type description of a model
         */
        explicit DetectorModel(std::string type) : type_(std::move(type)) {}
        /**
         * @brief Essential virtual destructor
         */
        virtual ~DetectorModel() = default;

        ///@{
        /**
         * @brief Use default copy and move behaviour
         */
        DetectorModel(const DetectorModel&) = default;
        DetectorModel& operator=(const DetectorModel&) = default;

        DetectorModel(DetectorModel&&) = default;
        DetectorModel& operator=(DetectorModel&&) = default;
        ///@}

        /**
         * @brief Get the type of the model
         * @return Model type
         */
        std::string getType() const { return type_; }
        /**
         * @brief Get size of the sensor
         * @return Size of the sensor
         */
        ROOT::Math::XYZVector getSensorSize() const { return sensor_size_; }

        /**
         * @brief Get size of the sensor in x-direction
         * @return Length in x
         */
        // TODO [doc] This method should be removed (or renamed)
        double getSensorSizeX() const { return sensor_size_.x(); };
        /**
         * @brief Get size of the sensor in x-direction
         * @return Length in y
         */
        // TODO [doc] This method should be removed (or renamed)
        double getSensorSizeY() const { return sensor_size_.y(); };
        /**
         * @brief Get size of the sensor in x-direction
         * @return Length in z
         */
        // TODO [doc] This method should be removed (or renamed)
        double getSensorSizeZ() const { return sensor_size_.z(); };

        /**
         * @brief Set the size of the sensor
         * @param val Size of the sensor
         */
        void setSensorSize(ROOT::Math::XYZVector val) { sensor_size_ = std::move(val); }

        /**
         * @brief Get number of pixel (replicated blocks in general sense)
         * @return List of two dimensional pixels
         */
        // NOTE: default to 1, no copied units
        virtual ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<int>> getNPixels() const {
            return ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<int>>(1, 1);
        }
        /**
         * @brief Get number of pixel in y-direction
         * @return Number of pixels in y-direction
         */
        virtual int getNPixelsX() const { return 1; };
        /**
         * @brief Get number of pixel in x-direction
         * @return Number of pixels in x-direction
         */
        virtual int getNPixelsY() const { return 1; };

        /**
         * @brief Get size of a single pixel
         * @return Size of a pixel
         */
        virtual ROOT::Math::XYVector getPixelSize() const {
            return ROOT::Math::XYVector(getSensorSizeX() / getNPixelsX(), getSensorSizeY() / getNPixelsY());
        }
        /**
         * @brief Get size of a single pixel in x-direction
         * @return Length in x
         */
        virtual double getPixelSizeX() const { return getPixelSize().x(); };
        /**
         * @brief Get size of a single pixel in y-direction
         * @return Length in y
         */
        virtual double getPixelSizeY() const { return getPixelSize().y(); };

        /**
         * @brief Get starting coordinate in x-direction of sensor in local frame of derived model
         * @return Minimum x-coordinate
         */
        virtual double getSensorMinX() const = 0;
        /**
         * @brief Get starting coordinate in y-direction of sensor in local frame of derived model
         * @return Minimum y-coordinate
         */
        virtual double getSensorMinY() const = 0;
        /**
         * @brief Get starting coordinate in z-direction of sensor in local frame of derived model
         * @return Minimum z-coordinate
         */
        virtual double getSensorMinZ() const = 0;

        /**
         * @brief Get coordinate of center of chip in local frame of derived model
         * @note It can be a bit counter intuitive that this is not always (0, 0, 0)
         *
         * The center coordinate corresponds to the \ref Detector::getPosition "position" in the global frame.
         */
        virtual ROOT::Math::XYZPoint getCenter() const = 0;

    private:
        std::string type_;

        ROOT::Math::XYZVector sensor_size_;
    };
} // namespace allpix

#endif // ALLPIX_GEOMETRY_MANAGER_H
